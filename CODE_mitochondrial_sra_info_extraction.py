#!/usr/bin/env python3
__version__ = 'b_thapamagar@mail.fhsu.edu|2025-03-10'

import argparse
import datetime
import json
import logging
import os
import bs4
import coloredlogs
import lxml
import pandas as pd
from CODE_mitochondrial_pubmed_article_extraction import PubmedInteract
from Bio import Entrez
import requests

def extract_library_info_from_sra_content(fulltext_etree: lxml.etree._Element):
    library_item = {}
    library_descriptor_list = fulltext_etree.xpath('//LIBRARY_DESCRIPTOR')

    if(len(library_descriptor_list) > 0):
    # Iterate through all child elements and print their tags and text
        for element in library_descriptor_list[0].iterchildren():
            key = element.tag
            value= ""
            if(element.tag == "LIBRARY_LAYOUT"):
                for layout_element in element.iterchildren():
                    value += layout_element.tag
            else:
                value = element.text
            library_item[key] = value
            
    platform_list = fulltext_etree.xpath('//PLATFORM')
    if(len(platform_list) > 0):
        platform = platform_list[0]
        instrument_model_list = platform.xpath('//INSTRUMENT_MODEL')
        if(len(instrument_model_list) > 0):
            instrument_model = instrument_model_list[0]
            library_item[instrument_model.tag] = instrument_model.text

    return library_item

def extract_SRA_info_from_ENA(ena_accession_id):
    # URL for POST request for European Nucleotide Archive
    post_url = 'https://www.ebi.ac.uk/ena/browser/api/xml'

    # Headers
    headers = {
        'accept': 'application/xml',
        'Content-Type': 'application/json'
    }

    # Data to be sent in the POST request
    request_body = {
        "accessions": [
            ena_accession_id
        ],
        "expanded": True,
        "annotationOnly": True,
        "lineLimit": 0,
        "download": False,
        "gzip": True,
        "set": True,
        "includeLinks": True,
        "range": "string",
        "complement": True
    }
    # Send the POST request
    response = requests.post(post_url, json=request_body, headers=headers)
    return response

def main(args):
    email = args.mail
    verbose = args.verbose
    
    formatted_datetime = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    if(not os.path.isdir("./log")):
        os.mkdir("./log")
    logger_filename = f"./log/{formatted_datetime}_mitochondrial_sra_sample_extraction.log"
    logging.basicConfig(filename=logger_filename, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO, logger=log)
    
    pubmed_interact = PubmedInteract(email= email, logger=log)
    
    # DATA_pubmed_info_with_bioproj_and_sra.csv file contains the information of the articles from PubMed with BioProject and SRA Code 
    # (Extracted by CODE_mitochondrial_pubmed_article_extraction.py and manually by Bryce)
    
    pubmed_info_with_bioproject = pd.read_csv("DATA_Combined_DATA_pubmed_info_with_bioproj_and_sra.csv")
    pubmed_info_with_bioproject = pubmed_info_with_bioproject.rename(columns={"European Nucleotide Archive": "ENA", "SRA Code": "SRA_Code"})
    print(pubmed_info_with_bioproject.info())
    
    pubmed_info_with_bioproject = pubmed_info_with_bioproject.fillna("")
    
    pubmed_dataframe = pd.DataFrame(columns=["Title", "BioProject_uid", "BioProject", "SRA_Code", "Pubmed_ID", "PMC_ID", "SRA_Count", "SRA_Id_list", "Is_in_Genbank"])
    
    for row in pubmed_info_with_bioproject.itertuples():
        log.info(row)
        pubmed_item = {
            "Title": row.TITLE,
            "BioProject": row.ENA,
            "SRA_Code": row.SRA_Code,
            "BioProject_uid": "",
            "Is_in_Genbank": False,
        }
        pubmed_id = ""
        pubmed_id= pubmed_interact.lookup_pubmed_id_by_title(row.TITLE)
        log.info(f"pubmed_id: {pubmed_id}")
        
        pmc_id = ""
        if(pubmed_id != ""):
            pmc_id = pubmed_interact.get_pmc_id_by_pubmed_id(pubmed_id)
            if(pmc_id != ""):
                log.info(f"pmc_id: {pmc_id}")
            else:
                log.error(f"pmc_id not found for pubmed_id {pubmed_id}")
        else:
            log.error(f"pubmed_id not found for title: {row.TITLE}")
        
        pubmed_item["Pubmed_ID"] = pubmed_id
        pubmed_item["PMC_ID"] = pmc_id
        
        if(row.ENA != ""):
            bio_proj_esearch_handle = Entrez.esearch(db="bioproject", term=row.ENA)
            bio_proj_esearch_result = Entrez.read(bio_proj_esearch_handle)
            if(len(bio_proj_esearch_result["IdList"]) > 0):
                pubmed_item["BioProject_uid"] = bio_proj_esearch_result["IdList"][0]
        elif(row.ENA == "" and row.SRA_Code != ""):
            sra_code = ""
            if("-" in row.SRA_Code):
                sra_records = row.SRA_Code.split("-")
                if(len(sra_records) > 0):
                    sra_code = sra_records[0]
            else:
                sra_code = row.SRA_Code
                
            if("SAMN" in sra_code):
                bio_sample_esearch_handle = Entrez.esearch(db="biosample", term=sra_code)
                bio_sample_result = Entrez.read(bio_sample_esearch_handle)
                bio_proj_elink_handle = Entrez.elink(dbfrom="biosample", db="bioproject", id=bio_sample_result["IdList"][0])
                bio_proj_elink_record = Entrez.read(bio_proj_elink_handle)
                pubmed_item["BioProject_uid"] = bio_proj_elink_record[0]["LinkSetDb"][0]["Link"][0]["Id"]
                
            else:
                sra_esearch_handle = Entrez.esearch(db="sra", term=sra_code)
                sra_search_result = Entrez.read(sra_esearch_handle)
                bio_proj_elink_handle = Entrez.elink(dbfrom="sra", db="bioproject", id=sra_search_result["IdList"][0])
                bio_proj_elink_record = Entrez.read(bio_proj_elink_handle)
                pubmed_item["BioProject_uid"] = bio_proj_elink_record[0]["LinkSetDb"][0]["Link"][0]["Id"]
            
            if(pubmed_item["BioProject_uid"] != ""):
                bio_proj_efetch_handle = Entrez.efetch(db="bioproject", id=pubmed_item["BioProject_uid"], retmode="xml")
                bioproj_etree = lxml.etree.fromstring(bio_proj_efetch_handle.read())
                archive_tag_list = bioproj_etree.xpath("//ArchiveID")
                if(len(archive_tag_list) > 0):
                    pubmed_item["BioProject"] = archive_tag_list[0].get("accession")
        if(pubmed_item["BioProject_uid"] != ""):    
            try:
                bio_proj_elink_handle = Entrez.elink(dbfrom="bioproject", db="sra", id=pubmed_item["BioProject_uid"])
                bio_proj_elink_record = Entrez.read(bio_proj_elink_handle)
                sra_id_list = []
                for link_set in bio_proj_elink_record:
                    for link in link_set["LinkSetDb"]:
                        for obj in link["Link"]:
                            sra_id_list.append(obj["Id"])
                sra_id_list = list(set(sra_id_list))
                sra_len = len(sra_id_list)
                sra_ids = ",".join(sra_id_list)
                pubmed_item["SRA_Id_list"] = sra_ids
                log.info(f"SRA ID count: {sra_len}")
                pubmed_item["SRA_Count"] = sra_len
                pubmed_item["Is_in_Genbank"] = True
                # print("Count: ", len(bio_proj_elink_record))
                pubmed_dataframe.loc[len(pubmed_dataframe)] = pubmed_item
            except Exception as ex:
                log.error(f"Exception ecnountered: {ex}")
                continue
        else:
            if(row.ENA != ""):
                url = f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_experiment&includeAccessions={row.ENA}&format=json"
                # Send the GET request
                response = requests.get(url)
                
                # Check the status code of the response
                if response.status_code == 200:
                    response_json = response.json()
                    pubmed_item["SRA_Count"] = len(response_json)
                    sra_id_list = []
                    for record in response_json:
                        sra_id_list.append(record["experiment_accession"])
                    sra_ids = ",".join(sra_id_list)
                    pubmed_item["SRA_Id_list"] = sra_ids
                    log.info("Bioproject found in European Nucleotide Archive.")
                else:
                    log.error("BioProject accession number not found.")
                    
            pubmed_dataframe.loc[len(pubmed_dataframe)] = pubmed_item  
    pubmed_dataframe.to_csv("DATA_pubmed_info_with_bioproj_sra_sample_new.csv", index=False)
    log.info("SRA Id Extraction Process completed.")
    
    pubmed_dataframe = pubmed_dataframe.fillna("")
    if(len(pubmed_dataframe) <= 0):
        pubmed_dataframe = pd.read_csv("DATA_pubmed_info_with_bioproj_sra_sample_new.csv",
                                    dtype={
                                        'BioProject_uid': 'string',
                                        'BioProject': 'string',
                                        'SRA_Id_list': 'string',
                                        'Is_in_Genbank': 'bool'}).fillna("")
    
    data = []
    for record in pubmed_dataframe.itertuples():
        sra_id_list =  []
        sra_id_list = record.SRA_Id_list.split(",")
        fulltext_soup = bs4.BeautifulSoup("", "html.parser")
        
        if(record.PMC_ID != "" and len(sra_id_list) > 0):
            article_content = pubmed_interact.get_complete_article_by_pmc_id(record.PMC_ID)
            fulltext_soup = bs4.BeautifulSoup(article_content, 'html.parser') 
            
        for sra_id in sra_id_list:
            log.info(f"Processing {record.BioProject}: {sra_id}")
            sra_item = {
                "BioProject": record.BioProject,
                "BioProject_uid": record.BioProject_uid,
                "SRA_Id": sra_id,
                "Sample_Title": ""
            }
            
            if(record.Is_in_Genbank):
                sra_fetch_handle = Entrez.efetch(db="sra", id=sra_id)
                sra_content = sra_fetch_handle.read()
                fulltext_etree = lxml.etree.fromstring(sra_content)
                title_tag_list = fulltext_etree.xpath("//TITLE")
                if(len(title_tag_list) > 0):
                    sra_item["Title"] = title_tag_list[0].text
                library_item = extract_library_info_from_sra_content(fulltext_etree)
                sra_item["Library"] = library_item
                
                # Extracting Sample Title from Pool Tag
                pool_tag_list = fulltext_etree.xpath("//Pool")
                if(len(pool_tag_list) > 0):
                    member_tag = pool_tag_list[0].find("Member")
                    sample_title = member_tag.get("sample_title")
                    sra_item["Original_Sample_Title"] = sample_title
                    if(sample_title):
                        if("_" in sample_title):
                            sample_title_splitted = sample_title.split("_")
                            if(len(sample_title_splitted) > 0):
                                sra_item["Sample_Title"] = sample_title_splitted[0]
                        else:
                            sra_item["Sample_Title"] = sample_title
                if(sra_item["Sample_Title"] != "" and fulltext_soup.contents):
                    log.info(f"Sample Title: {sra_item['Sample_Title']}")
                    target_content = None
                    target_content_list = fulltext_soup.find_all(string= lambda text: text and sra_item["Sample_Title"].lower() in text.lower())
                    if(len(target_content_list) > 0):
                        for target_content_item in target_content_list:
                            target_contentitem_tag = target_content_item.find_parent()
                            if((target_contentitem_tag.name == "td") or (target_contentitem_tag.name == "p")):
                                target_content = target_content_item
                                break
                    if(target_content):
                        target_tag = target_content.find_parent()
                        if(target_tag.name == "td"):
                            header_tag = target_tag.find_parent().find_parent().find_parent().find("thead")
                            table_column_name_list = []
                            for table_header_tag in header_tag.find_all("th"):
                                table_column_name_list.append(table_header_tag.text)
                            
                            if(len(table_column_name_list) == 0):
                                for table_header_tag in header_tag.find_all("td"):
                                    table_column_name_list.append(table_header_tag.text)
                            
                            data_row_tag = target_tag.find_parent()
                            data_row_list = []
                            for data_tag in data_row_tag.find_all("td"):
                                data_row_list.append(data_tag.text)
                            
                            table_obj = {}
                            for i in range(len(table_column_name_list)):
                                table_obj[table_column_name_list[i]] = data_row_list[i]    
                            
                            sra_item["Matched_Information"] = table_obj
                            
                            log.info(f"Matched Information: for {sra_item['BioProject']} - {sra_item['SRA_Id']} in tablular format.")
                        elif(target_tag.name == "p"):
                            sra_item["Matched_Information"] = target_tag.text
                            
                            log.info(f"Matched Information: for {sra_item['BioProject']} - {sra_item['SRA_Id']} in paragraph.")
                        else:
                            log.info(f"{target_tag.name} Case yet to handled but Matched Information: for {sra_item['BioProject']} - {sra_item['SRA_Id']}.")
                    else:
                        log.info("Sample Title not found in full text.")
                else:
                    if(not fulltext_soup.contents):
                        log.info("Full text not found.")
                    if(sra_item["Sample_Title"] == ""):
                        log.info("Sample Title not present.")
            else:
                response = extract_SRA_info_from_ENA(sra_id)
                if(response.status_code == 200):
                    fulltext_etree = lxml.etree.fromstring(response.text)
                    library_item = extract_library_info_from_sra_content(fulltext_etree)
                    sra_item["Library"] = library_item
                    log.info("Library information is extracted.")
                log.info("BioProject not in Genbank and ENA does not contain Sample Name.")
            data.append(sra_item)
        # else:
        #     log.info(f"PMC_ID: {record.PMC_ID} or SRA_ID: {record.SRA_Id_list} not found for {record.BioProject}")
    output_file_name = "DATA_sra_info_in_pmc.json"
    with open(output_file_name, "w") as file:
        json.dump(data, file, indent=4)
    log.info(f"Data written to {output_file_name}")
    log.info("Finished.")


if(__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Author|Version: '+__version__)
    parser.add_argument("--mail", "-m", type=str, required=True,
                        help="Your email address (needed for querying NCBI PubMed via Entrez)")
    
    parser.add_argument("--verbose", "-v", action="store_true", required=False, 
                        default=True, help="(Optional) Enable verbose logging")
    args = parser.parse_args()
    main(args)