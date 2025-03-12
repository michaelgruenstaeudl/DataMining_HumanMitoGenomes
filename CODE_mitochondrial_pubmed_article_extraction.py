#!/usr/bin/env python3
__version__ = 'b_thapamagar@mail.fhsu.edu|2025-02-18'

import argparse
from collections import namedtuple
import os
import Bio
import logging
import Bio.Entrez
import coloredlogs
import urllib
import lxml
import bs4
import json
from pathlib import Path
import pandas as pd
import re
import datetime

substrings: list = [
        "NCBI SRA",
        "Sequence Read Archive",
        "www.ncbi.nlm.gov/sra",
        "NCBI Sequence Read Archive",
        "European Nucleotide Archive",
        "ENA"
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
    "availability of data and material"
]

class LXMLops:

    def __init__(self, etree):
        self.etree = etree

    def remove_expendable(self):
        '''Remove unnecessary XML sections from full text'''
        xml_document = self.etree.find('.//document')
        # STEP 1. Removing figures, tables and any backmatter parts
        for passage in self.etree.findall('.//passage'):
            if passage.find('infon[@key="section_type"]').text in ['FIG', 'TABLE', 'COMP_INT', 'AUTH_CONT', 'ACK_FUND']:
                xml_document.remove(passage)  
        # STEP 2. Removing references
        for passage in self.etree.findall('.//passage'):
            if passage.find('infon[@key="type"]').text == 'ref':
                xml_document.remove(passage)

    def extract_all_text(self):
        '''Extract all text of the full text by paragraph'''
        paragraphs = []
        for passage in self.etree.findall('.//passage'):
            header = passage.find('infon[@key="section_type"]').text
            if header not in paragraphs:
                paragraphs.append(header)
            main_text = passage.find('text').text
            if main_text:
                paragraphs.append(main_text)
        #for paragr in paragraphs:
        #    for paragraph in re.findall(regex_for_sentence_delin, paragr):
        #        paragraphs.append(paragraph)
        return paragraphs

class PubmedInteract:
    def __init__(self, email, logger: logging.Logger):
        """
        Initializes the instance with the provided email and logger.

        Args:
            email (str): The email address to be used with NCBI Entrez.
            logger (logging.Logger): A logger instance for logging messages.
        """
        self.email = email
        Bio.Entrez.email = email
        self.logger = logger
        
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
    
        search_query = Bio.Entrez.esearch(db= "pubmed", 
                                        term= f"{title}[TITLE]",  
                                        sort='relevance')
        result = Bio.Entrez.read(search_query)
        
        if(int(result["Count"]) == 0):
            search_query = Bio.Entrez.esearch(db= "pubmed", 
                                            term= f"{title}", 
                                            sort='relevance')
            result = Bio.Entrez.read(search_query)
        
        if(int(result["Count"]) == 0):
            search_query = Bio.Entrez.esearch(db= "pubmed", 
                                            term= f"{title[: int(len(title)/2)]}",
                                            sort='relevance')
            result = Bio.Entrez.read(search_query)
        
        return result
    
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
        if(int(search_result["Count"]) > 0):
            for id in search_result["IdList"]:
                fetch_handle = Bio.Entrez.esummary(db= "pubmed", id= id)
                pubmed_result = Bio.Entrez.read(fetch_handle)
                if(title.lower() in pubmed_result[0]["Title"].lower()):
                    pubmed_id = pubmed_result[0]["Id"]
                    break
        
        if(pubmed_id == ""): 
            encoded_string = urllib.parse.quote(title)
            url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_string}"
            bioc_handle = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            pubmed_response = urllib.request.urlopen(bioc_handle).read()
            
            pubmed_soup = bs4.BeautifulSoup(pubmed_response, 'html.parser') 
            
            pubmed_citation_tag = (pubmed_soup.find("meta",attrs= {"name": "citation_pmid"}))
            
            if(pubmed_citation_tag):
                pubmed_id = pubmed_citation_tag["content"]
            
            if(pubmed_id == ""): 
                matching_citation_tag_list = pubmed_soup.find_all(name="section", attrs={"class" : "matching-citations search-results-list"})
                for docsum_tag in matching_citation_tag_list:
                    a_tag = docsum_tag.find(name="a", attrs={"class": "docsum-title"})
                    docsum_title = ''.join(a_tag.stripped_strings)
                    if(title in docsum_title):
                        try:
                            pubmed_id = a_tag["data-ga-label"]
                            break
                        except Exception as ex:
                            pubmed_id = ""
            
            if(pubmed_id == ""):
                displayed_uids_tag = (pubmed_soup.find("meta",attrs= {"name": "log_displayeduids"}))
                if(displayed_uids_tag):
                    pubmed_id_list = displayed_uids_tag["content"].split(",")
                    for id in pubmed_id_list:
                        fetch_handle = Bio.Entrez.esummary(db= "pubmed", id= id)
                        pubmed_result = Bio.Entrez.read(fetch_handle)
                        if(title.lower() in pubmed_result[0]["Title"].lower()):
                            pubmed_id = pubmed_result[0]["Id"]
                            break
                            
        return pubmed_id
    
    def fetch_pubmed_by_id(self, pubmed_id):
        pubmed_efetch_handle = Bio.Entrez.efetch(db="pubmed", id=pubmed_id)
        pubmed_result = Bio.Entrez.read(pubmed_efetch_handle)
        return pubmed_result
    
    def extract_url_to_full_article_by_id(self, pubmed_id):
        
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        pubmed_handle = urllib.request.Request(pubmed_url, headers={'User-Agent': 'Mozilla/5.0'})
        pubmed_article = urllib.request.urlopen(pubmed_handle).read()
        pubmed_soup = bs4.BeautifulSoup(pubmed_article, 'html.parser')
        
        full_text_link_div = pubmed_soup.find('div', class_='full-text-links-list')
        complete_article_link_list = []
        if(full_text_link_div):
            tag_a_list = (full_text_link_div.find_all('a'))
            for tag in tag_a_list:
                if(tag.has_attr('href')):
                    complete_article_link_list.append(tag["href"])
        return complete_article_link_list

    def get_pmc_id_by_pubmed_id(self, pubmed_id):
        '''Look up PubMedCentral ID from a PubMed ID'''
        handle = Bio.Entrez.elink(dbfrom='pubmed',
                                db='pmc',
                                linkname='pubmed_pmc',
                                id=pubmed_id,
                                retmode='text')
        
        result = Bio.Entrez.read(handle)
        try:
            pmcid = f"PMC{result[0]['LinkSetDb'][0]['Link'][0]['Id']}"
        except:
            pmcid = None
            
        return pmcid
        
    def get_complete_article_by_pmc_id(self, pmc_id):  
        '''
        Extract complete article from pubmed central based on PMC ID passed
        '''
        article_complete = None
        if(pmc_id != None):
            # article_url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmc_id}/unicode"
            # url_handle = urllib.request.urlopen(article_url)
            # article_complete = url_handle.read() 
            
            # if( "[Error] : No result can be found" in article_complete.decode('utf-8')):
            article_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc_id}/?report=reader"
            try:
                url_handle = urllib.request.Request(article_url, headers={'User-Agent': 'Mozilla/5.0'})
                article_complete = urllib.request.urlopen(url_handle).read()
            except:
                self.logger.warning(f"\t {pmc_id}: Retrieval unsuccessful")
                article_complete = None
                    
        else:
            self.logger.warning(f"No PMC Id available")
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
        all_paragraphs =[]
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
                        item = {
                            data_availability_item: all_paragraphs[i+1]
                        }
                        data_content_list.append(item)
        if '<!DOCTYPE html>' in str(article_complete):
            # Parse the HTML of the full text
            fulltext_soup = bs4.BeautifulSoup(article_complete, 'html.parser') 
            # cleaned_title = ''.join(char if char.isalnum() or char.isspace() else '' for char in article_title)
            # cleaned_title_from_html = ""
            # if(fulltext_soup.head.title.text):
            #     cleaned_title_from_html = ''.join(char if char.isalnum() or char.isspace() else '' for char in fulltext_soup.head.title.text)
            # if(cleaned_title.lower() in cleaned_title_from_html.lower()):
            article_content= fulltext_soup.find(name="section", attrs={"aria-label": "Article content"})
            for reflist_tag in article_content.find_all('section', class_ ="ref-list"):
                reflist_tag.decompose()
            
            target_tag_list = fulltext_soup.find_all(string= lambda text: text and text.lower() in data_availability_list)
            for target_tag in target_tag_list:
                parent_tag = target_tag.find_parent()
                if(parent_tag):
                    super_parent_tag = parent_tag.find_parent()
                    if(super_parent_tag):
                        data_content = super_parent_tag.find("p").text
                        item = {
                            target_tag.text: data_content
                        }
                        data_content_list.append(item)
            
            for paragraph in article_content.find_all("p"):
                text = ''.join(paragraph.stripped_strings)
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
                if(sub == "ENA"):
                    # start_idx = lower_para.find(f"\b{lower_sub}")
                    pattern = rf'\b{re.escape(lower_sub)}\b'
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

                    matching_paragraph_list.append({"paragraph":i+1,
                                    "substring": sub, 
                                    "content": content})
        
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
            
        if pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"] != []:
            pubmed_information["Published_Year"] = pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"]
            
        if "DataBankList" in pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]:
            pubmed_information["DataBankList"] = json.dumps(pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]["DataBankList"])
        
        complete_article_link_list = self.extract_url_to_full_article_by_id(pubmed_id)
        
        pubmed_information["Full_Article_URL"] = (", ".join(complete_article_link_list))
        
        pubmed_information["is_PMC"] = False
        for url_link in complete_article_link_list:
            if("https://pmc.ncbi.nlm.nih.gov/articles/pmid" in url_link):
                pubmed_information["is_PMC"] = True
                break
        return pubmed_information

def main(args):
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
    email = args.mail
    verbose = args.verbose
    file_path = args.filepath
    
    ### STEP 1. Set up logger
    # Configure the logging
    formatted_datetime = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    if(not os.path.isdir("./log")):
        os.mkdir("./log")
    logger_filename = f"./log/{formatted_datetime}_mitochondrial_pubmed_article_extraction.log"
    logging.basicConfig(filename=logger_filename, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO, logger=log)
    
    #STEP 2. Check if the file exists
    if Path(file_path).exists():
        log.info("File exists.")
    else:
        log.info("File does not exist.")
        return 
    
    pubmed_interact = PubmedInteract(email= email, logger= log)
    
    data_frame = pd.read_csv(file_path)
    title_list = data_frame["TITLE"].dropna().unique()
    
    pubmed_metadata = pd.DataFrame(columns= ["Pubmed_ID", "Title", "Published_Year", "DataBankList", "Full_Article_URL", "is_PMC", "Error"])
    
    #STEP3: Extraction pubmed metadata from title
    # for title in ["Neolithic phylogenetic continuity inferred from complete mitochondrial DNA sequences in a tribal population of Southern India"]:
    for title in title_list:
        item = {
            "Title": title
        }
        try:
            pubmed_id : str = ""
            log.info(f"querying PubMed for {title}")
            pubmed_id = pubmed_interact.lookup_pubmed_id_by_title(title)
            if(pubmed_id == ""):
                item["Error"] = "No pubmed id available."
                item["is_PMC"] = False
                pubmed_metadata.loc[len(pubmed_metadata)] = item
                log.info("No pubmed id available.")   
                continue  
            
            log.info(f"pubmed_id: {pubmed_id}")
            item["Pubmed_ID"] = pubmed_id
            additional_pubmed_info = pubmed_interact.get_pubmed_informations_by_pubmed_id(pubmed_id)
            item.update(additional_pubmed_info)
        except Exception as ex: 
            if("Error" not in item):
                item["Error"] = ""
            item["Error"] = f"\n {ex}" 
            log.info(f"Exception occured: {ex}")
            continue
        
        pubmed_metadata.loc[len(pubmed_metadata)] = item
    
    pubmed_metadata.to_csv("DATA_pubmed_metadata.csv", header=True, index=False)
    log.info("Pubmed ID extraction completed")
    
    ### STEP 4. Extracting full text and matching paragraphs
    if(len(pubmed_metadata)  == 0):
        pubmed_metadata = pd.read_csv("DATA_pubmed_metadata.csv",
                                    dtype={
                                        'Pubmed_ID': 'string',
                                        'DataBankList': 'string'})
    pubmed_metadata = pubmed_metadata.fillna("")
    
    log.info("Pubmed article extraction begins")
    data: list = []
    for row in pubmed_metadata[pubmed_metadata["Pubmed_ID"] != ""].itertuples():
        article_complete = None
        try:
            record = {
                "title": row.Title,
                "Pubmed_ID": row.Pubmed_ID,
                "DataBankList" : json.loads(row.DataBankList) if (row.DataBankList != "") else "",
                "URL":  row.Full_Article_URL,
                "Year": row.Published_Year
            }            
            pmc_id = pubmed_interact.get_pmc_id_by_pubmed_id(row.Pubmed_ID)
            if(pmc_id != None):
                log.info(f"pmd_id: {row.Pubmed_ID} and pmc_id: {pmc_id} \nretrieving complete article from PubMedCentral")
                record["pmc_id"] = pmc_id
                article_complete = pubmed_interact.get_complete_article_by_pmc_id(pmc_id)
            else:
                log.info(f"No PMC Id available for pubmed id {row.Pubmed_ID}")
                article_complete = None
            
            if(article_complete):
                try:
                    (all_paragraphs, data_content_list) = pubmed_interact.get_paragraph_list_from_pubmed_article(article_complete, row.Title)
                    record["DataContent"] = data_content_list
                    record["FullTextParagraph"] = all_paragraphs
                    matching_paragraph_list = pubmed_interact.get_matching_paragraphs_for_substrings(all_paragraphs, substrings)
                    if(matching_paragraph_list != []):
                        log.info(f"matched paragraphs for {record["title"]}")
                    record["MatchedParagraphs"] = matching_paragraph_list           
                except Exception as ex:
                    record["Error"] = f"Error encountered for {pmc_id} \n {ex}"
                    log.critical(f"Error encountered for {pmc_id} \n {ex}")
            else:
                record["Error"] = f"No content available"
            data.append(record)        
        except Exception as e:
            log.critical(f"Exception occured: {e}")
            
    matched_output_dict = [json_obj for json_obj in data if (json_obj.get('MatchedParagraphs') != None and json_obj.get('MatchedParagraphs') != [])]
    matched_output_file_name = "DATA_pubmed_records_with_data_source_info.json"
    with open(matched_output_file_name, "w") as file:
        json.dump(matched_output_dict, file, indent=4)
    output_file_name = "DATA_pubmed_records.json"
    with open(output_file_name, "w") as file:
        json.dump(data, file, indent=4)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Author|Version: '+__version__)
    parser.add_argument("--mail", "-m", type=str, required=True,
                        help="Your email address (needed for querying NCBI PubMed via Entrez)")
    
    parser.add_argument("--verbose", "-v", action="store_true", required=False, 
                        default=True, help="(Optional) Enable verbose logging")
    parser.add_argument("--filepath", "-f",  type=str, required=True, 
                        help="(Required) Filename which contain title")
    args = parser.parse_args()
    main(args= args)