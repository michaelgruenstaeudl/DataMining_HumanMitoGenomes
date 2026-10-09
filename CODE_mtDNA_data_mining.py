__version__ = "b_thapamagar@mail.fhsu.edu|2026-05-10"

import argparse
import json
import logging
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from difflib import SequenceMatcher
from http.client import IncompleteRead
from pathlib import Path
from urllib.error import HTTPError, URLError

import bs4
import coloredlogs
import lxml
import pandas as pd
from Bio import Entrez, SeqIO

# region Constant================================#
batch_size = 100  # Number of records to fetch in each batch
MAX_RETRIES = 5  # Maximum number of retries for fetching records
PUBMED_REQUEST_INTERVAL = 1.0  # Minimum seconds between paced PubMed requests
_last_request_time = None
substrings: list = [
    "NCBI SRA",
    "Sequence Read Archive",
    "www.ncbi.nlm.gov/sra",
    "NCBI Sequence Read Archive",
    "European Nucleotide Archive",
    "ENA",
]

data_availability_list: list = [
    "electronic-database information",
    "electronic database information",
    "associated data",
    "accession numbers",
    "data access",
    "data accessibility",
    "data availability",
    "availability of data",
    "data and code availability",
    "data availability statement",
    "availability of data and material",
]

repository_map = {
    "NCBI SRA": "NCBI Sequence Read Archive",
    "Sequence Read Archive": "NCBI Sequence Read Archive",
    "www.ncbi.nlm.gov/sra": "NCBI Sequence Read Archive",
    "NCBI Sequence Read Archive": "NCBI Sequence Read Archive",
    "European Nucleotide Archive": "European Nucleotide Archive",
    "ENA": "European Nucleotide Archive",
}

accession_patterns = {
    "NCBI Sequence Read Archive": r"\b(?:SRA|SRP|SRX|SRR|SRS|ERP|ERX|ERR|ERS)\d+",
    "European Nucleotide Archive": r"\b(?:PRJNA|PRJEA|PRJEB|ERP|ERX|ERR|ERS)\d+",
    "NCBI BioSample": r"\b(?:SAMN|SAMEA)\d+(?:[-\u2013](?:SAMN|SAMEA)\d+)?\b",
}
# endregion


# region Helper methods
def _retry_request(
    func, max_retries, logger=None, operation=None, request_interval=None
):
    global _last_request_time

    logger = logger or logging.getLogger(__name__)
    operation = operation or "network request"
    method, separator, operation_name = operation.partition(": ")
    context = f"[{method}] [{operation_name}]" if separator else f"[{operation}]"

    for attempt in range(1, max_retries + 1):
        try:
            if request_interval is not None and _last_request_time is not None:
                elapsed = time.monotonic() - _last_request_time
                pacing_delay = request_interval - elapsed
                if pacing_delay > 0:
                    logger.info(
                        f"{context} Waiting {pacing_delay:.2f} seconds "
                        "before PubMed request"
                    )
                    time.sleep(pacing_delay)

            if request_interval is not None and request_interval > 0:
                _last_request_time = time.monotonic()
            return func()
        except HTTPError as e:
            logger.error(
                f"{context} Network error on attempt " f"{attempt}/{max_retries}: {e}"
            )
            if attempt == max_retries:
                raise
            if e.code == 429:
                retry_after = e.headers.get("Retry-After") if e.headers else None

                if retry_after and retry_after.isdigit():
                    delay = max(int(retry_after), 90)
                else:
                    delay = 90 * attempt
            else:
                delay = 2 ** (attempt - 1)

            logger.info(f"{context} Waiting {delay} seconds before retrying")
            time.sleep(delay)
        except (IncompleteRead, URLError, OSError) as e:
            logger.error(
                f"{context} Network error on attempt " f"{attempt}/{max_retries}: {e}"
            )
            if attempt == max_retries:
                raise
            time.sleep(2 ** (attempt - 1))


def configure_logging(verbose):
    formatted_datetime = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M")
    log_directory = Path("./log")
    log_directory.mkdir(parents=True, exist_ok=True)
    logger_filename = (
        log_directory
        / f"{formatted_datetime}_UTC_mitochondrial_pubmed_article_extraction.log"
    )
    logging.basicConfig(
        filename=logger_filename,
        level=logging.DEBUG,
        format="%(asctime)s UTC - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)
    coloredlogs.install(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.DEBUG if verbose else logging.INFO,
        logger=logger,
    )
    return logger


# endregion


# Methods for fetching nucleotide summary records in batches
def fetch_batch(start, end, nucleotide_esearch_output, logger):
    return _retry_request(
        lambda: _fetch_nucleotide_batch(start, end, nucleotide_esearch_output),
        max_retries=MAX_RETRIES,
        logger=logger,
        operation="fetch_batch: Entrez.efetch nucleotide GenBank batch",
    )


def _fetch_nucleotide_batch(start, end, nucleotide_esearch_output):
    handle = Entrez.efetch(
        db="nucleotide",
        rettype="gb",
        retmode="text",
        id=",".join(nucleotide_esearch_output["IdList"][start:end]),
    )
    batch_records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    return batch_records


def extract_nucleotide_detailed_metadata_information(output_directory, logger):
    directory = Path(output_directory)
    nucleotide_metadata_output_file_path = (
        directory / "DATA_Nucleotide_detailed_metadata_records.csv"
    )
    search_term = """Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] 
                    NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"""

    nucleotide_esearch_output = _retry_request(
        lambda: Entrez.read(
            Entrez.esearch(
                db="Nucleotide",
                term=search_term,
                usehistory="y",
                retmax=100000,
            )
        ),
        max_retries=MAX_RETRIES,
        logger=logger,
        operation=(
            "extract_nucleotide_detailed_metadata_information: "
            "Entrez.esearch nucleotide records"
        ),
    )

    for start in range(0, len(nucleotide_esearch_output["IdList"]), batch_size):
        data_records = pd.DataFrame(
            columns=[
                "UID",
                "AccessionID",
                "BioProject",
                "BioSample",
                "VERSION",
                "ORGANISM",
                "SEQ_LEN",
                "CREATE_DATE",
                "AUTHORS",
                "TITLE",
                "REFERENCE",
                "TAXONOMY",
            ]
        )

        try:
            print(
                f"Fetching records in batch {start // batch_size + 1}: {batch_size} records"
            )

            end = min(start + batch_size, len(nucleotide_esearch_output["IdList"]))
            batch_records = fetch_batch(start, end, nucleotide_esearch_output, logger)

            print("Fetching completed and data extraction starts:")
            for record in batch_records:

                item = {
                    # "UID": uid,
                    "AccessionID": record.name,
                    "VERSION": record.annotations["sequence_version"],
                    "ORGANISM": record.annotations["organism"],
                    "SEQ_LEN": len(record.seq),
                    "DATE": record.annotations["date"],
                    "AUTHORS": record.annotations["references"][0].authors,
                    "TITLE": record.annotations["references"][0].title,
                    "REFERENCE": record.annotations["references"][0].journal,
                    "TAXONOMY": ";".join(record.annotations["taxonomy"]),
                }

                if record.dbxrefs != []:
                    # parsed_data = []
                    for string in record.dbxrefs:
                        key, value = string.split(":")
                        if key == "BioProject" or key == "BioSample":
                            item[key] = value
                        # parsed_data.append({key: value})

                data_records.loc[len(data_records)] = item

            if start == 0:
                data_records.to_csv(
                    nucleotide_metadata_output_file_path,
                    mode="w",
                    header=True,
                    index=False,
                )
            else:
                data_records.to_csv(
                    nucleotide_metadata_output_file_path,
                    mode="a",
                    header=False,
                    index=False,
                )

            print("Records saved successfully")
            print(
                f"Fetched batch {start // batch_size + 1}: {len(batch_records)} records"
            )

        except (KeyError, IndexError, ValueError, OSError) as e:
            data_records.to_csv(nucleotide_metadata_output_file_path, index=False)
            print("Records saved because of exception encountered")
            print(f"An error occurred on nucleotide metadata information: {e}")
    return data_records


def fetch_nucleotide_summary_batch(start, end, nucleotide_esearch_output, logger):
    batch = nucleotide_esearch_output["IdList"][start:end]
    batch_str = ",".join(batch)

    # Fetch summary data for the batch
    return _retry_request(
        lambda: Entrez.read(
            Entrez.esummary(db="nucleotide", id=batch_str, retmode="xml")
        ),
        max_retries=MAX_RETRIES,
        logger=logger,
        operation=(
            "fetch_nucleotide_summary_batch: " "Entrez.esummary nucleotide records"
        ),
    )


def extract_nucleotide_metadata_information(output_directory, logger):
    directory = Path(output_directory)
    nucleotide_metadata_output_file_path = (
        directory / "DATA_Nucleotide_Summary_records.csv"
    )

    search_term = """Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] 
                NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"""

    nucleotide_esearch_output = _retry_request(
        lambda: Entrez.read(
            Entrez.esearch(
                db="Nucleotide",
                term=search_term,
                usehistory="y",
                retmax=100000,
            )
        ),
        max_retries=MAX_RETRIES,
        logger=logger,
        operation=(
            "extract_nucleotide_metadata_information: "
            "Entrez.esearch nucleotide records"
        ),
    )
    summary_data_list = []
    try:
        for start in range(0, len(nucleotide_esearch_output["IdList"]), batch_size):
            end = min(start + batch_size, len(nucleotide_esearch_output["IdList"]))
            batch_records = fetch_nucleotide_summary_batch(
                start, end, nucleotide_esearch_output, logger
            )
            summary_data_list.extend(batch_records)  # Append batch_records to records
            print(
                f"Fetched batch {start // batch_size + 1}: {len(batch_records)} records"
            )

        # Convert records to DataFrame
        df = pd.DataFrame(summary_data_list)
        df.to_csv(nucleotide_metadata_output_file_path, index=False)
        print("Nucleotide summary records saved successfully")
    except (KeyError, OSError, ValueError) as e:
        print(f"An error occurred in nucleotide metadata extraction: {e}")

    return df


# Methods for fetching SRA records metadata in batches


def fetch_sra_metadata_batch(start, end, sra_esearch_output, logger):
    batch = sra_esearch_output["IdList"][start:end]
    batch_str = ",".join(batch)

    # Fetch summary data for the batch
    return _retry_request(
        lambda: Entrez.read(Entrez.esummary(db="sra", id=batch_str, retmode="xml")),
        max_retries=MAX_RETRIES,
        logger=logger,
        operation="fetch_sra_metadata_batch: Entrez.esummary SRA records",
    )


def extract_sra_metadata_batch(output_directory, logger):
    directory = Path(output_directory)
    sra_metadata_output_file_path = directory / "DATA_SRA_Summary_records.csv"
    search_term = """(human[organism] OR \"homo sapiens\"[organism]) AND (\"mitochondrial\"[title] or mitochondrion[TITLE])"""
    sra_esearch_output = _retry_request(
        lambda: Entrez.read(
            Entrez.esearch(
                db="sra",
                term=search_term,
                usehistory="y",
                retmax=100000,
            )
        ),
        max_retries=MAX_RETRIES,
        logger=logger,
        operation="extract_sra_metadata_batch: Entrez.esearch SRA records",
    )
    summary_data_list = []
    try:
        for start in range(0, len(sra_esearch_output["IdList"]), batch_size):
            end = min(start + batch_size, len(sra_esearch_output["IdList"]))
            batch_records = fetch_sra_metadata_batch(
                start, end, sra_esearch_output, logger
            )
            summary_data_list.extend(batch_records)  # Append batch_records to records
            print(
                f"Fetched batch {start // batch_size + 1}: {len(batch_records)} records"
            )

        # Convert records to DataFrame
        df = pd.DataFrame(summary_data_list)
        df.to_csv(sra_metadata_output_file_path, index=False)
        print("SRA summary records saved successfully")
    except (KeyError, OSError, ValueError) as e:
        print(f"An error occurred in SRA metadata extraction: {e}")


# Classes for pubmed data mining


class LXMLops:

    def __init__(self, etree):
        self.etree = etree

    def remove_expendable(self):
        """Remove unnecessary XML sections from full text"""
        xml_document = self.etree.find(".//document")
        # STEP 1. Removing figures, tables and any backmatter parts
        for passage in self.etree.findall(".//passage"):
            if passage.find('infon[@key="section_type"]').text in [
                "FIG",
                "TABLE",
                "COMP_INT",
                "AUTH_CONT",
                "ACK_FUND",
            ]:
                xml_document.remove(passage)
        # STEP 2. Removing references
        for passage in self.etree.findall(".//passage"):
            if passage.find('infon[@key="type"]').text == "ref":
                xml_document.remove(passage)

    def extract_all_text(self):
        """Extract all text of the full text by paragraph"""
        paragraphs = []
        for passage in self.etree.findall(".//passage"):
            header = passage.find('infon[@key="section_type"]').text
            if header not in paragraphs:
                paragraphs.append(header)
            main_text = passage.find("text").text
            if main_text:
                paragraphs.append(main_text)
        return paragraphs


class PubmedInteract:
    def __init__(
        self,
        email,
        logger: logging.Logger,
        request_interval=PUBMED_REQUEST_INTERVAL,
    ):
        """
        Initializes the instance with the provided email and logger.

        Args:
            email (str): The email address to be used with NCBI Entrez.
            logger (logging.Logger): A logger instance for logging messages.
        """
        self.email = email
        self.logger = logger
        self._request_interval = request_interval

    def _retry_request(self, func, operation, request_interval=None):
        return _retry_request(
            func,
            max_retries=MAX_RETRIES,
            logger=self.logger,
            operation=f"{self.__class__.__name__}.{operation}",
            request_interval=request_interval,
        )

    def build_title_query(self, title):
        # Extract words while ignoring punctuation
        words = re.findall(r"\b[\w]+\b", title)

        # PubMed stop words that don't help identify the article
        stop_words = {
            "a",
            "an",
            "and",
            "are",
            "as",
            "at",
            "be",
            "by",
            "for",
            "from",
            "in",
            "is",
            "of",
            "on",
            "or",
            "the",
            "to",
            "with",
        }

        words = [word for word in words if word.lower() not in stop_words]

        return " AND ".join(f'"{word}"[Title]' for word in words)

    def _normalize_title(self, title):
        title = unicodedata.normalize("NFKC", title)
        title = title.lower()

        # Treat different quotation marks as equivalent by removing punctuation.
        title = re.sub(r"[^\w\s]", " ", title)

        # Collapse multiple spaces.
        title = re.sub(r"\s+", " ", title).strip()

        return title

    def _calculate_title_similarity(self, title1, title2):
        """
        Calculate normalized similarity between two article titles.

        Returns a value between 0.0 and 1.0.
        """
        # normalized_title1 = self._normalize_title(title1)
        # normalized_title2 = self._normalize_title(title2)

        return SequenceMatcher(None, title1, title2).ratio()

    def search_pubmed_by_title(self, title):
        """
        Search PubMed for articles by title.
        This method searches the PubMed database for articles that match the given title.
        It performs three levels of search:
        1. Exact match for the title.
        2. General search for the title.
        3. Search for the first half of the title.
        Parameters:
        title (str): The title of the article to search for.
        Returns:
        dict: A dictionary containing the search results from PubMed.
        """

        def search():
            term = self.build_title_query(title)
            result = self._retry_request(
                lambda: Entrez.read(
                    Entrez.esearch(db="pubmed", term=term, sort="relevance")
                ),
                operation="search_pubmed_by_title: Entrez.esearch title-term search",
            )

            if int(result["Count"]) == 0:
                term = self._normalize_title(title)
                result = self._retry_request(
                    lambda: Entrez.read(
                        Entrez.esearch(db="pubmed", term=f"{term}", sort="relevance")
                    ),
                    operation="search_pubmed_by_title: Entrez.esearch normalized title search",
                )

            if int(result["Count"]) == 0:
                result = self._retry_request(
                    lambda: Entrez.read(
                        Entrez.esearch(
                            db="pubmed",
                            term=f"{title[: int(len(title)/2)]}",
                            sort="relevance",
                        )
                    ),
                    operation=(
                        "search_pubmed_by_title: "
                        "Entrez.esearch first-half title search"
                    ),
                )

            return result

        return search()

    def lookup_pubmed_id_by_title(self, title):
        """
        Lookup the PubMed ID for a given article title.
        This method attempts to find the PubMed ID associated with a given article title by first searching the PubMed database using the Entrez API. If no results are found, it then performs a web scraping operation on the PubMed website to locate the PubMed ID.
        Args:
            title (str): The title of the article to search for.
        Returns:
            str: The PubMed ID of the article if found, otherwise an empty string.
        """
        pubmed_id: str = ""
        search_result = self.search_pubmed_by_title(title)
        if int(search_result["Count"]) > 0:
            exact_match_found = False
            best_match_id = None
            best_similarity = 0.0
            for id in search_result["IdList"]:
                pubmed_result = self._retry_request(
                    lambda id=id: Entrez.read(Entrez.esummary(db="pubmed", id=id)),
                    operation="lookup_pubmed_id_by_title: Entrez.esummary title match",
                )

                normalized_title = self._normalize_title(title)
                normalized_pubmed_title = self._normalize_title(
                    pubmed_result[0]["Title"]
                )

                if normalized_title == normalized_pubmed_title:
                    pubmed_id = pubmed_result[0]["Id"]
                    exact_match_found = True
                    break

                similarity = self._calculate_title_similarity(
                    normalized_title, normalized_pubmed_title
                )

                if similarity > best_similarity:
                    best_similarity = similarity
                    best_match_id = pubmed_result[0]["Id"]

            # 3. Use similarity only if no exact match was found
            if not exact_match_found and best_similarity >= 0.90:
                pubmed_id = best_match_id

        if pubmed_id == "":
            encoded_string = urllib.parse.quote(title)
            url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_string}"
            pubmed_response = self._retry_request(
                lambda: urllib.request.urlopen(
                    urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                ).read(),
                operation=(
                    "lookup_pubmed_id_by_title: " "urllib.request.urlopen search page"
                ),
                request_interval=self._request_interval,
            )
            pubmed_soup = bs4.BeautifulSoup(pubmed_response, "html.parser")

            pubmed_citation_tag = pubmed_soup.find(
                "meta", attrs={"name": "citation_pmid"}
            )
            if pubmed_citation_tag:
                pubmed_id = pubmed_citation_tag["content"]

            if pubmed_id == "":
                matching_citation_tag_list = pubmed_soup.find_all(
                    name="section",
                    attrs={"class": "matching-citations search-results-list"},
                )
                for docsum_tag in matching_citation_tag_list:
                    a_tag = docsum_tag.find(name="a", attrs={"class": "docsum-title"})
                    docsum_title = "".join(a_tag.stripped_strings)
                    if title in docsum_title and a_tag is not None:
                        pubmed_id = a_tag.get("data-ga-label", "")
                        if pubmed_id:
                            break

            if pubmed_id == "":
                displayed_uids_tag = pubmed_soup.find(
                    "meta", attrs={"name": "log_displayeduids"}
                )
                if displayed_uids_tag:
                    pubmed_id_list = displayed_uids_tag["content"].split(",")

                    for id in pubmed_id_list:
                        pubmed_result = self._retry_request(
                            lambda id=id: Entrez.read(
                                Entrez.esummary(db="pubmed", id=id)
                            ),
                            operation=(
                                "lookup_pubmed_id_by_title: "
                                "Entrez.esummary displayed UID title match"
                            ),
                        )
                        if title.lower() in pubmed_result[0]["Title"].lower():
                            pubmed_id = pubmed_result[0]["Id"]
                            break
        return pubmed_id

    def fetch_pubmed_by_id(self, pubmed_id):
        return self._retry_request(
            lambda: Entrez.read(Entrez.efetch(db="pubmed", id=pubmed_id)),
            operation="fetch_pubmed_by_id: Entrez.efetch PubMed record",
        )

    def extract_url_to_full_article_by_id(self, pubmed_id):

        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        pubmed_article = self._retry_request(
            lambda: urllib.request.urlopen(
                urllib.request.Request(
                    pubmed_url, headers={"User-Agent": "Mozilla/5.0"}
                )
            ).read(),
            operation=(
                "extract_url_to_full_article_by_id: "
                "urllib.request.urlopen PubMed page"
            ),
            request_interval=self._request_interval,
        )
        pubmed_soup = bs4.BeautifulSoup(pubmed_article, "html.parser")

        full_text_link_div = pubmed_soup.find("div", class_="full-text-links-list")
        complete_article_link_list = []
        if full_text_link_div:
            tag_a_list = full_text_link_div.find_all("a")
            for tag in tag_a_list:
                if tag.has_attr("href"):
                    complete_article_link_list.append(tag["href"])
        return complete_article_link_list

    def get_pmc_id_by_pubmed_id(self, pubmed_id):
        """Look up PubMedCentral ID from a PubMed ID"""
        result = self._retry_request(
            lambda: Entrez.read(
                Entrez.elink(
                    dbfrom="pubmed",
                    db="pmc",
                    linkname="pubmed_pmc",
                    id=pubmed_id,
                    retmode="text",
                )
            ),
            operation="get_pmc_id_by_pubmed_id: Entrez.elink PubMed-to-PMC lookup",
        )
        try:
            pmcid = f"PMC{result[0]['LinkSetDb'][0]['Link'][0]['Id']}"
        except (IndexError, KeyError):
            pmcid = None

        return pmcid

    def get_complete_article_by_pmc_id(self, pmc_id):
        """
        Extract complete article from pubmed central based on PMC ID passed
        """
        if pmc_id != None:
            # article_url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmc_id}/unicode"
            # url_handle = urllib.request.urlopen(article_url)
            # article_complete = url_handle.read()

            # if( "[Error] : No result can be found" in article_complete.decode('utf-8')):
            article_url = (
                f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc_id}/?report=reader"
            )
            try:
                article_complete = self._retry_request(
                    lambda: urllib.request.urlopen(
                        urllib.request.Request(
                            article_url, headers={"User-Agent": "Mozilla/5.0"}
                        )
                    ).read(),
                    operation=(
                        "get_complete_article_by_pmc_id: "
                        "urllib.request.urlopen PMC article"
                    ),
                    request_interval=self._request_interval,
                )
            except (HTTPError, URLError, OSError):
                self.logger.warning(f"\t {pmc_id}: Retrieval unsuccessful")
                return None

        else:
            self.logger.warning("No PMC Id available")
            article_complete = None

        return article_complete

    def get_paragraph_list_from_pubmed_article(self, article_complete, article_title):
        """
        Extracts paragraphs and specific data content from a PubMed article.
        This method processes the provided PubMed article, which can be in XML or HTML format,
        and extracts all paragraphs and specific data content based on predefined data availability items.
        Args:
            article_complete (str): The complete content of the PubMed article in XML or HTML format.
            article_title (str): The title of the PubMed article.
        Returns:
            tuple: A tuple containing:
                - all_paragraphs (list): A list of all paragraphs extracted from the article.
                - data_content_list (list): A list of dictionaries containing specific data content
                    extracted based on predefined data availability items.
        """
        all_paragraphs = []
        data_content_list = []
        if '<?xml version="1.0"' in str(article_complete):
            # Parse the XML of the full text
            fulltext_etree = lxml.etree.fromstring(article_complete)

            # Remove unnecessary sections from full text
            LXMLops(fulltext_etree).remove_expendable()
            # Extract all text of the full text by paragraph
            all_paragraphs = LXMLops(fulltext_etree).extract_all_text()
            for data_availability_item in data_availability_list:
                for i in range(len(all_paragraphs) - 1):
                    if all_paragraphs[i].lower() == data_availability_item:
                        item = {data_availability_item: all_paragraphs[i + 1]}
                        data_content_list.append(item)
        if "<!DOCTYPE html>" in str(article_complete):
            # Parse the HTML of the full text
            fulltext_soup = bs4.BeautifulSoup(article_complete, "html.parser")
            # cleaned_title = ''.join(char if char.isalnum() or char.isspace() else '' for char in article_title)
            # cleaned_title_from_html = ""
            # if(fulltext_soup.head.title.text):
            #     cleaned_title_from_html = ''.join(char if char.isalnum() or char.isspace() else '' for char in fulltext_soup.head.title.text)
            # if(cleaned_title.lower() in cleaned_title_from_html.lower()):
            article_content = fulltext_soup.find(
                name="section", attrs={"aria-label": "Article content"}
            )
            for reflist_tag in article_content.find_all("section", class_="ref-list"):
                reflist_tag.decompose()

            target_tag_list = fulltext_soup.find_all(
                string=lambda text: text and text.lower() in data_availability_list
            )
            for target_tag in target_tag_list:
                parent_tag = target_tag.find_parent()
                if parent_tag:
                    super_parent_tag = parent_tag.find_parent()
                    if super_parent_tag:
                        data_content = super_parent_tag.find("p").text
                        item = {target_tag.text: data_content}
                        data_content_list.append(item)

            for paragraph in article_content.find_all("p"):
                text = "".join(paragraph.stripped_strings)
                all_paragraphs.append(text)
            # else:
            #     self.logger.info("Different document pulled")
        return (all_paragraphs, data_content_list)

    def get_matching_paragraphs_for_substrings(self, paragraph_list, substring_list):
        """
        Find and extract paragraphs containing specified substrings.
        This method searches through a list of paragraphs and identifies those that contain any of the specified substrings.
        For each match, it extracts a portion of the paragraph surrounding the substring and stores the result.
        Args:
            paragraph_list (list of str): A list of paragraphs to search through.
            substring_list (list of str): A list of substrings to search for within the paragraphs.
        Returns:
            list of dict: A list of dictionaries, each containing:
                - "paragraph" (int): The index of the paragraph (1-based).
                - "substring" (str): The substring that was found.
                - "content" (str): A portion of the paragraph surrounding the found substring, with up to 100 characters before and after the substring.
        """
        matching_paragraph_list = []  # Store matching results
        for i, paragraph in enumerate(paragraph_list):
            for sub in substring_list:
                lower_para = paragraph.lower()
                lower_sub = sub.lower()
                if sub == "ENA":
                    # start_idx = lower_para.find(f"\b{lower_sub}")
                    pattern = rf"\b{re.escape(lower_sub)}\b"
                    # Find the starting index of the first match
                    match = re.search(pattern, lower_para)
                    start_idx = match.start() if match else -1

                else:
                    start_idx = lower_para.find(lower_sub)

                if start_idx != -1:  # If the substring is found
                    # Extract 100 chars before and after, ensuring we don't go out of bounds
                    start = max(0, start_idx - 100)
                    end = min(len(paragraph), start_idx + len(sub) + 100)
                    content = paragraph[start:end]

                    matching_paragraph_list.append(
                        {"paragraph": i + 1, "substring": sub, "content": content}
                    )

        return matching_paragraph_list

    def get_pubmed_informations_by_pubmed_id(self, pubmed_id):
        """
        Retrieve PubMed information for a given PubMed ID.
        Args:
            pubmed_id (str): The PubMed ID of the article to retrieve information for.
        Returns:
            dict: A dictionary containing the following keys:
                - "Published_Year" (str): The year the article was published.
                - "DataBankList" (str): A JSON string of the DataBankList associated with the article.
                - "Full_Article_URL" (str): A comma-separated string of URLs to the full article.
                - "is_PMC" (bool): A flag indicating whether the article is available in PMC (PubMed Central).
        """
        pubmed_information = {}
        pubmed_result = self.fetch_pubmed_by_id(pubmed_id)

        if (
            pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "ArticleDate"
            ]
            != []
        ):
            pubmed_information["Published_Year"] = pubmed_result["PubmedArticle"][0][
                "MedlineCitation"
            ]["Article"]["ArticleDate"][0]["Year"]

        if (
            "DataBankList"
            in pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]
        ):
            pubmed_information["DataBankList"] = json.dumps(
                pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"][
                    "DataBankList"
                ]
            )

        # this code always lead to http 429 (Too Many Requests) error. So, no url is extracted from pubmed instead we use elink to check pmcID availability
        # ==============================#
        # complete_article_link_list = self.extract_url_to_full_article_by_id(pubmed_id)

        # pubmed_information["Full_Article_URL"] = ", ".join(complete_article_link_list)

        # pubmed_information["is_PMC"] = False
        # for url_link in complete_article_link_list:
        #     if "https://pmc.ncbi.nlm.nih.gov/articles/pmid" in url_link:
        #         pubmed_information["is_PMC"] = True
        #         break
        # ==============================#
        # For now full_article_url is set as empty string
        pubmed_information["Full_Article_URL"] = ""
        pubmed_information["is_PMC"] = False

        return pubmed_information


# Pubmed Article Information Extraction


def extract_pubmed_article_information_by_title(
    args,
    output_directory: Path,
    title_list: list,
    pubmed_interact: PubmedInteract,
    logger,
):
    """
    Main function to extract PubMed article metadata and full text based on titles from a CSV file.
    Args:
        args (argparse.Namespace): Command-line arguments containing:
            - mail (str): Email address for PubMed API.
            - verbose (bool): Verbosity flag for logging.
            - filepath (str): Path to the input CSV file containing article titles.
    Steps:
        1. Set up logger for logging information and errors.
        2. Check if the input file exists.
        3. Extract PubMed metadata for each title in the CSV file.
        4. Save the extracted metadata to a CSV file.
        5. Extract full text and matching paragraphs from PubMed articles.
        6. Save the matched paragraphs and full data to JSON files.
    Returns:
        None
    """
    # email = args.mail
    # STEP 2. Check if the file exists
    # if Path(file_path).exists():
    #     log.info("File exists.")
    # else:
    #     log.info("File does not exist.")
    #     return
    # ------------------------------#

    # data_frame = pd.read_csv(file_path)

    pubmed_metadata = pd.DataFrame(
        columns=[
            "Pubmed_ID",
            "Title",
            "Published_Year",
            "DataBankList",
            "Full_Article_URL",
            "is_PMC",
            "Error",
        ]
    )

    # STEP3: Extraction pubmed metadata from title
    # for title in ["Neolithic phylogenetic continuity inferred from complete mitochondrial DNA sequences in a tribal population of Southern India"]:
    for i, title in enumerate(title_list):
        print(f"Processing title {i + 1} of {len(title_list)}: {title}")
        item = {"Title": title}
        try:
            pubmed_id: str = ""
            logger.info(f"querying PubMed for {title}")
            pubmed_id = pubmed_interact.lookup_pubmed_id_by_title(title)
            if pubmed_id == "":
                item["Error"] = "No pubmed id available."
                item["is_PMC"] = False
                pubmed_metadata.loc[len(pubmed_metadata)] = item
                logger.info("No pubmed id available.")
                continue

            logger.info(f"pubmed_id: {pubmed_id}")
            item["Pubmed_ID"] = pubmed_id
            additional_pubmed_info = (
                pubmed_interact.get_pubmed_informations_by_pubmed_id(pubmed_id)
            )
            item.update(additional_pubmed_info)
        except Exception as ex:
            if "Error" not in item:
                item["Error"] = ""
            item["Error"] = f"\n {ex}"
            logger.exception(
                "[extract_pubmed_article_information_by_title] "
                "Unexpected error while processing article"
            )
            continue

        pubmed_metadata.loc[len(pubmed_metadata)] = item

    pubmed_metadata.to_csv(
        output_directory / "DATA_pubmed_metadata.csv", header=True, index=False
    )
    logger.info("Pubmed ID extraction completed")

    ### STEP 4. Extracting full text and matching paragraphs
    if len(pubmed_metadata) == 0:
        pubmed_metadata = pd.read_csv(
            output_directory / "DATA_pubmed_metadata.csv",
            dtype={"Pubmed_ID": "string", "DataBankList": "string"},
        )
    pubmed_metadata = pubmed_metadata.fillna("")

    logger.info("Pubmed article extraction begins")
    data: list = []
    for row in pubmed_metadata[pubmed_metadata["Pubmed_ID"] != ""].itertuples():
        article_complete = None
        try:
            record = {
                "title": row.Title,
                "Pubmed_ID": row.Pubmed_ID,
                "DataBankList": (
                    json.loads(row.DataBankList) if (row.DataBankList != "") else ""
                ),
                "URL": row.Full_Article_URL,
                "Year": row.Published_Year,
            }
            pmc_id = pubmed_interact.get_pmc_id_by_pubmed_id(row.Pubmed_ID)
            if pmc_id != None:
                logger.info(
                    f"pmd_id: {row.Pubmed_ID} and pmc_id: {pmc_id} \nretrieving complete article from PubMedCentral"
                )
                record["pmc_id"] = pmc_id
                article_complete = pubmed_interact.get_complete_article_by_pmc_id(
                    pmc_id
                )
            else:
                logger.info(f"No PMC Id available for pubmed id {row.Pubmed_ID}")
                article_complete = None

            if article_complete:
                try:
                    all_paragraphs, data_content_list = (
                        pubmed_interact.get_paragraph_list_from_pubmed_article(
                            article_complete, row.Title
                        )
                    )
                    record["DataContent"] = data_content_list
                    record["FullTextParagraph"] = all_paragraphs
                    matching_paragraph_list = (
                        pubmed_interact.get_matching_paragraphs_for_substrings(
                            all_paragraphs, substrings
                        )
                    )
                    if matching_paragraph_list != []:
                        logger.info(f"matched paragraphs for {record['title']}")
                    record["MatchedParagraphs"] = matching_paragraph_list
                except Exception as ex:
                    record["Error"] = f"Error encountered for {pmc_id} \n {ex}"
                    # logger.critical(f"Error encountered for {pmc_id} \n {ex}")
                    logger.exception(
                        f"[extract_pubmed_full_text_by_pmc_id] for {pmc_id}"
                        "Unexpected error while processing article"
                    )
            else:
                record["Error"] = "No content available"
            data.append(record)
        except Exception:
            logger.exception(
                f"[pubmed_article_mining] for {pmc_id}"
                "Unexpected error while processing article"
            )

    matched_output_dict = [
        json_obj
        for json_obj in data
        if (
            json_obj.get("MatchedParagraphs") != None
            and json_obj.get("MatchedParagraphs") != []
        )
    ]
    matched_output_file_path = (
        output_directory / "DATA_pubmed_records_with_data_source_info.json"
    )
    with open(matched_output_file_path, "w") as file:
        json.dump(matched_output_dict, file, indent=4)
    output_file_path = output_directory / "DATA_pubmed_records.json"
    with open(output_file_path, "w") as file:
        json.dump(data, file, indent=4)

    return matched_output_dict


# Extract BioProject ID, BioSample ID, SRA ID, and ENA ID from the extracted pubmed article data content.


def extract_data_source_info_from_pubmed_article_data_content(
    output_directory: Path,
    matched_pubmed_record,
    logger,
):
    logger.info("Extracting data source information from pubmed article data content")

    results = []
    # for item in [d for d in matched_records if d.get("title") == "DNA analysis of an early modern human from Tianyuan Cave, China"]:
    for item in matched_pubmed_record:
        logger.info(f"Processing item: {item.get('title')}")

        logger.info("Processing data content")
        for data_content in item.get("DataContent", []):
            key, value = next(iter(data_content.items()))
            for repository, pattern in accession_patterns.items():
                accessions = re.findall(
                    pattern,
                    value,
                    flags=re.IGNORECASE,
                )

                results.append(
                    {
                        "title": item.get("title"),
                        "repository": repository,
                        "accessions": ", ".join(dict.fromkeys(accessions)),
                        "substring": key,
                    }
                )

        logger.info("Processing matched paragraphs")
        for matched_paragraph_item in item["MatchedParagraphs"]:
            substring = matched_paragraph_item.get("substring", "")
            content = matched_paragraph_item.get("content", "")
            repository = repository_map.get(substring)

            if repository is None:
                continue

            accessions = re.findall(
                accession_patterns[repository],
                content,
                flags=re.IGNORECASE,
            )

            results.append(
                {
                    "title": item.get("title"),
                    "repository": repository,
                    "accessions": ", ".join(dict.fromkeys(accessions)),
                    "substring": substring,
                }
            )

    results_df = pd.DataFrame(results)

    logger.info(
        "Filtering and prioritizing data source information based on accessions and repository"
    )
    # Step 1: Create a copy of the original DataFrame
    df = results_df.copy()

    # Step 2: Identify rows with non-empty accessions

    priority = {
        "Data Availability Statement": 1,
        "Data availability": 1,
        "Associated Data": 2,
        "ENA": 3,
    }

    # Step 3: Create a new column "_priority" based on the "substring" column
    df["_priority"] = df["substring"].map(priority).fillna(99)

    df["_has_accessions"] = df["accessions"].fillna("").str.strip().ne("")

    # Step 4: Identify rows from the European Nucleotide Archive
    df["_is_ena"] = df["repository"].eq("European Nucleotide Archive")

    # Step 5: Sort by title and selection priority
    df = df.sort_values(
        by=["title", "_has_accessions", "_priority", "_is_ena"],
        ascending=[True, False, True, False],
        kind="stable",
    )

    # Step 6: Keep the first row for each title
    df = df.drop_duplicates(
        subset=["title"],
        keep="first",
    )

    # Step 7: Remove temporary columns
    df = df.drop(columns=["_has_accessions", "_is_ena", "_priority"]).reset_index(
        drop=True
    )

    logger.info("Saving filtered data to CSV")
    df.to_csv(
        f"{output_directory}/DATA_pubmed_records_with_data_source_info_filtered.csv",
        index=False,
    )

    return df


# Mapping of Nucleotide to SRA records:
# def mapping():


def main(args):

    # region Portal for the entire data mining pipeline testing

    # Replaces step 1: Testing the data mining pipeline with pre-extracted nucleotide metadata and detailed metadata information
    # ------------------------------#
    # nucleotide_detailed_metadata_info = pd.read_csv(
    #     "test_output/DATA_Nucleotide_detailed_metadata_records.csv"
    # )
    # nucleotide_metadata_info = pd.read_csv(
    #     "test_output/DATA_Nucleotide_Summary_records.csv"
    # )
    # ------------------------------#

    # Replaces step 3: Testing the data mining pipeline with pre-extracted pubmed article information
    # ------------------------------#

    # matched_records = []
    # with open(
    #     "test_output/DATA_pubmed_records_with_data_source_info.json", "r"
    # ) as file:
    #     matched_records = json.load(file)
    # ------------------------------#

    # endregion

    ##############################################################
    ######## Data mining pipeline execution starts here ##########
    ##############################################################

    logger = configure_logging(args.verbose)
    directory = Path(args.output_directory)

    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)

    API_KEY = (
        args.api_key
    )  # Please extract api key from NCBI website (sign in and get access to API KEY)
    Entrez.email = (
        args.mail
    )  # Mention the email address same as which is used to sign in NCBI.
    Entrez.api_key = API_KEY

    # Step 1: Fetching nucleotide summary records for Homo sapiens complete mitochondrial genome sequences
    nucleotide_metadata_info = extract_nucleotide_metadata_information(
        args.output_directory, logger
    )
    nucleotide_detailed_metadata_info = (
        extract_nucleotide_detailed_metadata_information(args.output_directory, logger)
    )

    # Step 2: SRA records metadata extraction
    extract_sra_metadata_batch(args.output_directory)

    # Step 3: Pubmed article mining Mining

    title_list = nucleotide_detailed_metadata_info["TITLE"].dropna().unique()

    pubmed_interact = PubmedInteract(email=args.mail, logger=logger)
    matched_output_dict = extract_pubmed_article_information_by_title(
        args, directory, title_list, pubmed_interact, logger
    )

    # Step 4: Extracting BioProject ID, BioSample ID, SRA ID, and ENA ID from the extracted pubmed article data content.
    # This step would involve processing the matched_output_dict to extract the required IDs.
    extracted_data_source_info_df = (
        extract_data_source_info_from_pubmed_article_data_content(
            directory, matched_output_dict, logger
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Author|Version: " + __version__)
    parser.add_argument(
        "--mail",
        "-m",
        type=str,
        required=True,
        help="Your email address (needed for querying NCBI PubMed via Entrez)",
    )

    parser.add_argument(
        "--output_directory",
        "-o",
        type=str,
        required=True,
        help="The directory where the output files will be saved",
    )

    parser.add_argument(
        "--api_key",
        "-k",
        type=str,
        required=True,
        help="Your API key (needed for querying NCBI PubMed via Entrez)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        required=False,
        default=True,
        help="(Optional) Enable verbose logging",
    )
    args = parser.parse_args()
    main(args)
