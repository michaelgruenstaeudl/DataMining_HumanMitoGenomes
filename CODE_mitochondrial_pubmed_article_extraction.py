#!/usr/bin/env python3
__version__ = 'b_thapamagar@mail.fhsu.edu|2025-02-18'

import argparse
import Bio
import logging
import Bio.Entrez
import coloredlogs
from litMiningPubmed_conductMining import HTMLops, LXMLops, PubMedInteract
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


def main(args):
    email = args.mail
    verbose = args.verbose
    file_path = args.filepath
    
    ### STEP 2. Set up logger
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
    
    
    data: list = []
    
    Bio.Entrez.email = email
    
    data_frame = pd.read_csv(file_path)
    title_list = data_frame["TITLE"].dropna()
    
    for title in ["Sephardic signature in haplogroup T mitochondrial DNA"]:
        try:
            record = {
                "title": title
            }
            
            log.info(f"querying PubMed for {title}")
            search_query = Bio.Entrez.esearch(db= "pubmed", term= f"{title}[TITLE]")
            content = Bio.Entrez.read(search_query)
            pubmed_id : str
            
            if(int(content["Count"]) > 0):
                for id in content["IdList"]:
                    fetch_handle = Bio.Entrez.esummary(db= "pubmed", id= id)
                    pubmed_result = Bio.Entrez.read(fetch_handle)
                    if(title in pubmed_result[0]["Title"]):
                        pubmed_id = pubmed_result[0]["Id"]
                        break
                    
            record["publed_id"] = pubmed_id
            
            elink_handle = Bio.Entrez.elink(dbfrom= "pubmed", id= pubmed_id, db="pmc")
            elink_result = Bio.Entrez.read(elink_handle)
            
            if(elink_result[0]["LinkSetDb"] == []):
                record["Error"] = "No Pubmed central record available"
                log.info(f"No pmc record for {title}")
                continue
            
            pmc_id = elink_result[0]["LinkSetDb"][0]["Link"][0]["Id"]
            
            log.info(f"pmd_id: {pubmed_id} and pmc_id: PCC{pmc_id} \nretrieving complete article from PubMedCentral")    
            
            record["pmc_id"] = pmc_id
            
            bioc_url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC{pmc_id}/unicode"
        
            bioc_handle = urllib.request.urlopen(bioc_url)
            article_complete = bioc_handle.read()
            
            if(article_complete):
                try:
                    if '<?xml version="1.0"' in str(article_complete):
                        # Parse the XML of the full text
                        fulltext_etree = lxml.etree.fromstring(article_complete)
                        date = fulltext_etree.find(".//infon[@key='year']").text
                        record["Year"] = date
                        if((not date) or (int(date) < 2009)):
                            log.info("Article published before 2009 and it doesnot contain SRA records.")
                            continue
                        
                        # Remove unnecessary sections from full text
                        LXMLops(fulltext_etree).remove_expendable()
                        # Extract all text of the full text by paragraph
                        all_paragraphs = LXMLops(fulltext_etree).extract_all_text()
                        record["FullTextParagraph"] = all_paragraphs
                        results = []  # Store matching results

                        for i, paragraph in enumerate(all_paragraphs):  
                            for sub in substrings:
                                lower_para = paragraph.lower()
                                lower_sub = sub.lower()
                                start_idx = lower_para.find(lower_sub)

                                if start_idx != -1:  # If the substring is found
                                    # Extract 100 chars before and after, ensuring we don't go out of bounds
                                    start = max(0, start_idx - 100)
                                    end = min(len(paragraph), start_idx + len(sub) + 100)
                                    content = paragraph[start:end]

                                    results.append({"paragraph":i+1,
                                                    "substring": sub, 
                                                    "content": content})
                        
                        record["MatchedParagraphs"] = results           
                        
                        print("Fetched")
                        
                except Exception as ex:
                    record["Error"] = f"Error encountered for {pmc_id} \n {ex}"
                    log.critical(f"Error encountered for {pmc_id} \n {ex}")
            else:
                record["Error"] = f"No content available"
            
            data.append(record)        
        except Exception as e:
            log.critical(f"Exception occured: {e}")
            
    output_file_name = "output.json"
    with open(output_file_name, "w") as file:
        json.dump(data, file, indent=4)
        
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