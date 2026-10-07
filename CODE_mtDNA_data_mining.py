__version__ = "b_thapamagar@mail.fhsu.edu|2026-05-10"

import argparse
import json
import logging
import os
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from http.client import IncompleteRead
from pathlib import Path
from urllib.error import HTTPError, URLError

import bs4
import coloredlogs
import lxml
import pandas as pd
from Bio import Entrez, SeqIO

# Constant================================#
batch_size = 100  # Number of records to fetch in each batch
MAX_RETRIES = 5  # Maximum number of retries for fetching records
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
# ===============================#


def _retry_request(func, logger, max_retries=MAX_RETRIES):
    for attempt in range(max_retries):
        try:
            return func()
        except (IncompleteRead, HTTPError, URLError, OSError) as e:
            logger.error(
                f"Network error on attempt {attempt + 1}/{max_retries}: {e}"
            )
            if attempt == max_retries - 1:
                raise
            time.sleep(2**attempt)


# Methods for fetching nucleotide summary records in batches
def fetch_batch(start, end, nucleotide_esearch_output):
    logger = logging.getLogger(__name__)
    return _retry_request(
        lambda: _fetch_nucleotide_batch(start, end, nucleotide_esearch_output),
        logger,
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


def extract_nucleotide_detailed_metadata_information(output_directory):
    directory = Path(output_directory)
    nucleotide_metadata_output_file_path = (
        directory / "DATA_Nucleotide_detailed_metadata_records.csv"
    )
    search_term = """Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] 
                    NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"""

    handle = Entrez.esearch(
        db="Nucleotide",
        term=search_term,
        usehistory="y",
        retmax=100000,  # Total count was: 62173 so, to access all UIDs, set retmax to 100000
    )

    nucleotide_esearch_output = Entrez.read(handle)

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
            batch_records = fetch_batch(start, end, nucleotide_esearch_output)

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


def fetch_nucleotide_summary_batch(start, end, nucleotide_esearch_output):
    batch = nucleotide_esearch_output["IdList"][start:end]
    batch_str = ",".join(batch)

    # Fetch summary data for the batch
    fetch_handle = Entrez.esummary(db="nucleotide", id=batch_str, retmode="xml")
    summary_data = Entrez.read(fetch_handle)
    fetch_handle.close()
    return summary_data


def extract_nucleotide_metadata_information(output_directory):
    directory = Path(output_directory)
    nucleotide_metadata_output_file_path = (
        directory / "DATA_Nucleotide_Summary_records.csv"
    )

    search_term = """Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] 
                NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"""

    handle = Entrez.esearch(
        db="Nucleotide",
        term=search_term,
        usehistory="y",
        retmax=100000,  # Total count was: 62173 so, to access all UIDs, set retmax to 100000
    )

    nucleotide_esearch_output = Entrez.read(handle)
    summary_data_list = []
    try:
        for start in range(0, len(nucleotide_esearch_output["IdList"]), batch_size):
            end = min(start + batch_size, len(nucleotide_esearch_output["IdList"]))
            batch_records = fetch_nucleotide_summary_batch(
                start, end, nucleotide_esearch_output
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


def fetch_sra_metadata_batch(start, end, sra_esearch_output):
    batch = sra_esearch_output["IdList"][start:end]
    batch_str = ",".join(batch)

    # Fetch summary data for the batch
    fetch_handle = Entrez.esummary(db="sra", id=batch_str, retmode="xml")
    summary_data = Entrez.read(fetch_handle)
    fetch_handle.close()
    return summary_data


def extract_sra_metadata_batch(output_directory):
    directory = Path(output_directory)
    sra_metadata_output_file_path = directory / "DATA_SRA_Summary_records.csv"
    search_term = """(human[organism] OR \"homo sapiens\"[organism]) AND (\"mitochondrial\"[title] or mitochondrion[TITLE])"""
    search_handle = Entrez.esearch(
        db="sra",
        term=search_term,
        usehistory="y",
        retmax=100000,  # Total count was: 62173 so, to access all UIDs, set retmax to 100000
    )
    sra_esearch_output = Entrez.read(search_handle)
    summary_data_list = []
    try:
        for start in range(0, len(sra_esearch_output["IdList"]), batch_size):
            end = min(start + batch_size, len(sra_esearch_output["IdList"]))
            batch_records = fetch_sra_metadata_batch(start, end, sra_esearch_output)
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
    MAX_RETRIES = MAX_RETRIES

    def __init__(self, email, logger: logging.Logger):
        """
        Initializes the instance with the provided email and logger.

        Args:
            email (str): The email address to be used with NCBI Entrez.
            logger (logging.Logger): A logger instance for logging messages.
        """
        self.email = email
        Entrez.email = email
        self.logger = logger

    def _retry_request(self, func):
        return _retry_request(func, self.logger, self.MAX_RETRIES)

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
            result = self._retry_request(
                lambda: Entrez.read(
                    Entrez.esearch(
                        db="pubmed", term=f"{title}[TITLE]", sort="relevance"
                    )
                )
            )

            if int(result["Count"]) == 0:
                result = self._retry_request(
                    lambda: Entrez.read(
                        Entrez.esearch(
                            db="pubmed", term=f"{title}", sort="relevance"
                        )
                    )
                )

            if int(result["Count"]) == 0:
                result = self._retry_request(
                    lambda: Entrez.read(
                        Entrez.esearch(
                            db="pubmed",
                            term=f"{title[: int(len(title)/2)]}",
                            sort="relevance",
                        )
                    )
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
            for id in search_result["IdList"]:
                pubmed_result = self._retry_request(
                    lambda id=id: Entrez.read(Entrez.esummary(db="pubmed", id=id))
                )
                if title.lower() in pubmed_result[0]["Title"].lower():
                    pubmed_id = pubmed_result[0]["Id"]
                    break

        if pubmed_id == "":
            encoded_string = urllib.parse.quote(title)
            url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_string}"
            pubmed_response = self._retry_request(
                lambda: urllib.request.urlopen(
                    urllib.request.Request(
                        url, headers={"User-Agent": "Mozilla/5.0"}
                    )
                ).read()
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
                    # Flag to indicate if a matching title was found
                    found = False

                    for id in pubmed_id_list:
                        pubmed_result = self._retry_request(
                            lambda id=id: Entrez.read(
                                Entrez.esummary(db="pubmed", id=id)
                            )
                        )
                        if title.lower() in pubmed_result[0]["Title"].lower():
                            found = True
                            pubmed_id = pubmed_result[0]["Id"]
                            break

                        if found:
                            break
        return pubmed_id

    def fetch_pubmed_by_id(self, pubmed_id):
        return self._retry_request(
            lambda: Entrez.read(Entrez.efetch(db="pubmed", id=pubmed_id))
        )

    def extract_url_to_full_article_by_id(self, pubmed_id):

        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        pubmed_article = self._retry_request(
            lambda: urllib.request.urlopen(
                urllib.request.Request(
                    pubmed_url, headers={"User-Agent": "Mozilla/5.0"}
                )
            ).read()
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
            )
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
                    ).read()
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

        complete_article_link_list = self.extract_url_to_full_article_by_id(pubmed_id)

        pubmed_information["Full_Article_URL"] = ", ".join(complete_article_link_list)

        pubmed_information["is_PMC"] = False
        for url_link in complete_article_link_list:
            if "https://pmc.ncbi.nlm.nih.gov/articles/pmid" in url_link:
                pubmed_information["is_PMC"] = True
                break
        return pubmed_information


# Pubmed Article Information Extraction


def extract_pubmed_article_information_by_title(
    args,
    output_directory: Path,
    nucleotide_metadata_df,
    pubmed_interact: PubmedInteract,
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
    verbose = args.verbose
    # file_path = args.filepath

    ### STEP 1. Set up logger
    # Configure the logging
    formatted_datetime = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M")
    if not os.path.isdir("./log"):
        os.mkdir("./log")
    logger_filename = (
        f"./log/{formatted_datetime}_UTC_mitochondrial_pubmed_article_extraction.log"
    )
    logging.basicConfig(
        filename=logger_filename,
        level=logging.DEBUG,
        format="%(asctime)s UTC - %(levelname)s - %(message)s",
    )
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.DEBUG,
            logger=log,
        )
    else:
        coloredlogs.install(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.INFO,
            logger=log,
        )

    # STEP 2. Check if the file exists
    # if Path(file_path).exists():
    #     log.info("File exists.")
    # else:
    #     log.info("File does not exist.")
    #     return

    # data_frame = pd.read_csv(file_path)
    title_list = nucleotide_metadata_df["TITLE"].dropna().unique()

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
            log.info(f"querying PubMed for {title}")
            pubmed_id = pubmed_interact.lookup_pubmed_id_by_title(title)
            if pubmed_id == "":
                item["Error"] = "No pubmed id available."
                item["is_PMC"] = False
                pubmed_metadata.loc[len(pubmed_metadata)] = item
                log.info("No pubmed id available.")
                continue

            log.info(f"pubmed_id: {pubmed_id}")
            item["Pubmed_ID"] = pubmed_id
            additional_pubmed_info = (
                pubmed_interact.get_pubmed_informations_by_pubmed_id(pubmed_id)
            )
            item.update(additional_pubmed_info)
        except Exception as ex:
            if "Error" not in item:
                item["Error"] = ""
            item["Error"] = f"\n {ex}"
            log.info(f"Exception occured: {ex}")
            continue

        pubmed_metadata.loc[len(pubmed_metadata)] = item

    pubmed_metadata.to_csv(
        output_directory / "DATA_pubmed_metadata.csv", header=True, index=False
    )
    log.info("Pubmed ID extraction completed")

    ### STEP 4. Extracting full text and matching paragraphs
    if len(pubmed_metadata) == 0:
        pubmed_metadata = pd.read_csv(
            output_directory / "DATA_pubmed_metadata.csv",
            dtype={"Pubmed_ID": "string", "DataBankList": "string"},
        )
    pubmed_metadata = pubmed_metadata.fillna("")

    log.info("Pubmed article extraction begins")
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
                log.info(
                    f"pmd_id: {row.Pubmed_ID} and pmc_id: {pmc_id} \nretrieving complete article from PubMedCentral"
                )
                record["pmc_id"] = pmc_id
                article_complete = pubmed_interact.get_complete_article_by_pmc_id(
                    pmc_id
                )
            else:
                log.info(f"No PMC Id available for pubmed id {row.Pubmed_ID}")
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
                        log.info(f"matched paragraphs for {record["title"]}")
                    record["MatchedParagraphs"] = matching_paragraph_list
                except Exception as ex:
                    record["Error"] = f"Error encountered for {pmc_id} \n {ex}"
                    log.critical(f"Error encountered for {pmc_id} \n {ex}")
            else:
                record["Error"] = "No content available"
            data.append(record)
        except Exception as e:
            log.critical(f"Exception occured: {e}")

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


# Mapping of Nucleotide to SRA records:
# def mapping():


def main(args):

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
    # nucleotide_metadata_info = extract_nucleotide_metadata_information(
    #     args.output_directory
    # )
    # nucleotide_detailed_metadata_info = (
    #     extract_nucleotide_detailed_metadata_information(args.output_directory)
    # )

    nucleotide_detailed_metadata_info = pd.read_csv(
        "test_output/DATA_Nucleotide_detailed_metadata_records.csv"
    )
    # nucleotide_metadata_info = pd.read_csv(
    #     "test_output/DATA_Nucleotide_Summary_records.csv"
    # )

    # Step 2: SRA records metadata extraction
    # extract_sra_metadata_batch(args.output_directory)

    # Step 3: Pubmed article mining Mining
    log = logging.getLogger(__name__)
    pubmed_interact = PubmedInteract(email=args.mail, logger=log)
    extract_pubmed_article_information_by_title(
        args, directory, nucleotide_detailed_metadata_info, pubmed_interact
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
