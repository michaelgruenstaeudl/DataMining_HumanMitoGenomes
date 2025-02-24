#!/usr/bin/env python3
__version__ = 'b_thapamagar@mail.fhsu.edu|2025-02-18'

import argparse
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


substrings: list = [
        "NCBI SRA",
        "Sequence Read Archive",
        "www.ncbi.nlm.gov/sra",
        "NCBI Sequence Read Archive",
        "European Nucleotide Archive",
    ]

class LXMLops:

    def __init__(self, etree):
        self.etree = etree

    def remove_expendable(self):
        '''Remove unnecessary XML sections from full text'''
        xml_document = self.etree.find('.//document')
        # STEP 1. Removing figures, tables and any backmatter parts
        for passage in self.etree.findall('.//passage'):
            if passage.find('infon[@key="section_type"]').text in ['FIG', 'TABLE', 'COMP_INT', 'AUTH_CONT', 'SUPPL', 'ACK_FUND']:
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
    def __init__(self, email):
        self.email = email
        Bio.Entrez.email = email
        
    def search_pubmed_by_title(self, title):
        '''Search PubMed via a query'''
        Bio.Entrez.email = self.email
    
        search_query = Bio.Entrez.esearch(db= "pubmed", term= f"{title}[TITLE]")
        result = Bio.Entrez.read(search_query)
        
        if(int(result["Count"]) == 0):
            search_query = Bio.Entrez.esearch(db= "pubmed", term= f"{title}")
            result = Bio.Entrez.read(search_query)
        
        if(int(result["Count"]) == 0):
            search_query = Bio.Entrez.esearch(db= "pubmed", term= f"{title[: int(len(title)/2)]}")
            result = Bio.Entrez.read(search_query)
        
        return result
    
    def lookup_pubmed_id_by_title(self, title):
        pubmed_id: str = ""
        search_result = self.search_pubmed_by_title(title)
        if(int(search_result["Count"]) > 0):
                for id in search_result["IdList"]:
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

def main(args):
    email = args.mail
    verbose = args.verbose
    file_path = args.filepath
    
    ### STEP 2. Set up logger
    # Configure the logging
    logging.basicConfig(filename='app.log', level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO, logger=log)
    
    # Check if the file exists
    if Path(file_path).exists():
        log.info("File exists.")
    else:
        log.info("File does not exist.")
        return 
    
    
    # Bio.Entrez.email = email
    pubmed_interact = PubmedInteract(email= email)
    
    data_frame = pd.read_csv(file_path)
    title_list = data_frame["TITLE"].dropna()
    
    pubmed_metadata = pd.DataFrame(columns= ["Pubmed_ID", "Title", "Published Year", "DataBankList", "Full Article URL", "is PMC", "Error"])
    
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
                item["is PMC"] = False
                pubmed_metadata.loc[len(pubmed_metadata)] = item
                log.info("No pubmed id available.")      
                continue  
            
            log.info(f"pubmed_id: {pubmed_id}")
            
            pubmed_result = pubmed_interact.fetch_pubmed_by_id(pubmed_id)
            
            item["Pubmed_ID"] = pubmed_id
            
            if pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"] != []:
                item["Published Year"] = pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"]
                
            if "DataBankList" in pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]:
                item["DataBankList"] = json.dumps(pubmed_result["PubmedArticle"][0]["MedlineCitation"]["Article"]["DataBankList"])
            
            complete_article_link_list = pubmed_interact.extract_url_to_full_article_by_id(pubmed_id)
            
            item["Full Article URL"] = (", ".join(complete_article_link_list))
            
            item["is PMC"] = False
            for url_link in complete_article_link_list:
                if("https://pmc.ncbi.nlm.nih.gov/articles/pmid" in url_link):
                    item["is PMC"] = True
                    break
            
        except Exception as ex: 
            if("Error" not in item):
                item["Error"] = ""
            item["Error"] = f"\n {ex}" 
            log.info(f"Exception occured: {ex}")
            continue
        
        pubmed_metadata.loc[len(pubmed_metadata)] = item
    
    pubmed_metadata.to_csv("pubmed_metadata.csv", header=True, index=False)
        
        
    #     try:
    #         record = {
    #             "title": title
    #         }
            
            
    #         record["publed_id"] = pubmed_id
            
    #         elink_handle = Bio.Entrez.elink(dbfrom= "pubmed", id= pubmed_id, db="pmc")
    #         elink_result = Bio.Entrez.read(elink_handle)
            
    #         if(elink_result[0]["LinkSetDb"] == []):
    #             record["Error"] = "No Pubmed central record available"
    #             log.info(f"No pmc record for {title}")
    #             continue
            
    #         pmc_id = elink_result[0]["LinkSetDb"][0]["Link"][0]["Id"]
            
    #         log.info(f"pmd_id: {pubmed_id} and pmc_id: PCC{pmc_id} \nretrieving complete article from PubMedCentral")    
            
    #         record["pmc_id"] = pmc_id
            
    #         bioc_url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC{pmc_id}/unicode"
        
    #         bioc_handle = urllib.request.urlopen(bioc_url)
    #         article_complete = bioc_handle.read()
            
    #         if(article_complete):
    #             try:
    #                 if '<?xml version="1.0"' in str(article_complete):
    #                     # Parse the XML of the full text
    #                     fulltext_etree = lxml.etree.fromstring(article_complete)
    #                     date = fulltext_etree.find(".//infon[@key='year']").text
    #                     record["Year"] = date
    #                     if((not date) or (int(date) < 2009)):
    #                         log.info("Article published before 2009 and it doesnot contain SRA records.")
    #                         continue
                        
    #                     # Remove unnecessary sections from full text
    #                     LXMLops(fulltext_etree).remove_expendable()
    #                     # Extract all text of the full text by paragraph
    #                     all_paragraphs = LXMLops(fulltext_etree).extract_all_text()
    #                     record["FullTextParagraph"] = all_paragraphs
    #                     results = []  # Store matching results

    #                     for i, paragraph in enumerate(all_paragraphs):  
    #                         for sub in substrings:
    #                             lower_para = paragraph.lower()
    #                             lower_sub = sub.lower()
    #                             start_idx = lower_para.find(lower_sub)

    #                             if start_idx != -1:  # If the substring is found
    #                                 # Extract 100 chars before and after, ensuring we don't go out of bounds
    #                                 start = max(0, start_idx - 100)
    #                                 end = min(len(paragraph), start_idx + len(sub) + 100)
    #                                 content = paragraph[start:end]

    #                                 results.append({"paragraph":i+1,
    #                                                 "substring": sub, 
    #                                                 "content": content})
                        
    #                     record["MatchedParagraphs"] = results           
                        
    #                     print("Fetched")
                        
    #             except Exception as ex:
    #                 record["Error"] = f"Error encountered for {pmc_id} \n {ex}"
    #                 log.critical(f"Error encountered for {pmc_id} \n {ex}")
    #         else:
    #             record["Error"] = f"No content available"
            
    #         data.append(record)        
    #     except Exception as e:
    #         log.critical(f"Exception occured: {e}")
            
    # output_file_name = "output.json"
    # with open(output_file_name, "w") as file:
    #     json.dump(data, file, indent=4)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Author|Version: '+__version__)
    parser.add_argument("--mail", "-m", type=str, required=True,
                        help="Your email address (needed for querying NCBI PubMed via Entrez)")
    
    parser.add_argument("--verbose", "-v", action="store_true", required=False, 
                        default=True, help="(Optional) Enable verbose logging")
    # parser.add_argument("--title", "-t", required=True, help="Paper title to query in PubMed")
    parser.add_argument("--filepath", "-f",  type=str, required=True, 
                        help="(Required) Filename which contain title")
    args = parser.parse_args()
    main(args= args)