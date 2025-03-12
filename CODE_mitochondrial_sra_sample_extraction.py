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

def main():
    formatted_datetime = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    if(not os.path.isdir("./log")):
        os.mkdir("./log")
    logger_filename = f"./log/{formatted_datetime}_mitochondrial_sra_sample_extraction.log"
    logging.basicConfig(filename=logger_filename, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger(__name__)
    # if verbose:
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    
    pubmed_interact = PubmedInteract(email= "b_thapamagar@mail.fhsu.edu", logger=log)
    pubmed_info_with_bioproject = pd.read_csv("DATA_pubmed_info_with_bioproj_and_sra.csv")
    # print(pubmed_info_with_bioproject.head())
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
            
            
    pubmed_dataframe.to_csv("DATA_pubmed_info_with_bioproj_sra_sample.csv", index=False)
    log.info("SRA Id Extraction Process completed.")
    
    pubmed_dataframe = pd.read_csv("DATA_pubmed_info_with_bioproj_sra_sample.csv",
                                    dtype={
                                        'BioProject_uid': 'string',
                                        'BioProject': 'string',
                                        'SRA_Id_list': 'string',
                                        'Is_in_Genbank': 'bool'}).fillna("")
    if(len(pubmed_dataframe) <= 0):
        return
    
    data = []
    for record in pubmed_dataframe.itertuples():
        sra_id_list =  []
        sra_id_list = record.SRA_Id_list.split(",")
        if(record.PMC_ID != "" and len(sra_id_list) > 0):
            article_content = pubmed_interact.get_complete_article_by_pmc_id(record.PMC_ID)
            fulltext_soup = bs4.BeautifulSoup(article_content, 'html.parser') 
            
            for sra_id in sra_id_list:
                log.info(f"Processing {record.BioProject}: {sra_id}")
                sra_item = {
                    "BioProject": record.BioProject,
                    "SRA_Id": sra_id,
                    "Sample_Title": ""
                }
                
                if(record.Is_in_Genbank):
                    sra_fetch_handle = Entrez.efetch(db="sra", id=sra_id)
                    sra_content = sra_fetch_handle.read()
                    fulltext_etree = lxml.etree.fromstring(sra_content)

                    pool_tag_list = fulltext_etree.xpath("//Pool")
                    if(len(pool_tag_list) > 0):
                        member_tag = pool_tag_list[0].find("Member")
                        sample_title = member_tag.get("sample_title")
                        if(sample_title):
                            sra_item["Sample_Title"] = sample_title
                    if(sra_item["Sample_Title"] != ""):
                        log.info(f"Sample Title: {sra_item['Sample_Title']}")
                        target_content = fulltext_soup.find(string= lambda text: text and sra_item["Sample_Title"].lower() in text.lower())
                        if(target_content):
                            target_tag = target_content.find_parent()
                            if(target_tag.name == "td"):
                                header_tag = target_tag.find_parent().find_parent().find_parent().find("thead")
                                table_column_name_list = []
                                for table_header_tag in header_tag.find_all("th"):
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
                        log.info("Sample Title not present.")
                else:
                    log.info("BioProject not in Genbank.")
                data.append(sra_item)
        else:
            log.info(f"PMC_ID: {record.PMC_ID} or SRA_ID: {record.SRA_Id_list} not found for {record.BioProject}")
    output_file_name = "DATA_sra_info_in_pmc.json"
    with open(output_file_name, "w") as file:
        json.dump(data, file, indent=4)
    log.info(f"Data written to {output_file_name}")
    log.info("Finished.")


if(__name__=="__main__"):
    main()