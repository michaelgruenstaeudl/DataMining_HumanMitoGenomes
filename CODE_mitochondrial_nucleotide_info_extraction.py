#!/usr/bin/env python3
__version__ = 'b_thapamagar@mail.fhsu.edu|2025-02-18'

import argparse
import datetime
import json
import logging
import os
import coloredlogs
from Bio import Entrez, SeqIO
import pandas as pd

class NucleotideInteract():
    def __init__(self, email, logger):
        self.email = email
        Entrez.email = self.email
        self.logger = logger
        
    def fetch_nucleotide_info(self, genbank_accession_id):
        self.logger.info(f"Fetching nucleotide information for {genbank_accession_id}")
        handle = Entrez.efetch(db="nucleotide", id=genbank_accession_id, rettype="gb", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record
    
    def extract_source_info_from_nucleotide_info(self, nucleotide_info):
        source_info = {}
        source_obj = next(filter(lambda x: x["GBFeature_key"] == "source", nucleotide_info[0]["GBSeq_feature-table"]), None)
        if source_obj:
            for qualifier in source_obj["GBFeature_quals"]:
                source_info[qualifier["GBQualifier_name"]] = qualifier["GBQualifier_value"]
        return source_info

def main(args):
    formatted_datetime = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    
    if(not os.path.isdir("./log/Nucleotide")):
        if(not os.path.isdir("./log")):
            os.mkdir("./log")
        os.mkdir("./log/Nucleotide")
    logger_filename = f"./log/Nucleotide/{formatted_datetime}_mitochondrial_nucleotide_info_extraction.log"
    logging.basicConfig(filename=logger_filename, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    log = logging.getLogger(__name__)
    # if verbose:
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    
    nucleotide_interact = NucleotideInteract(email="b_thapamagar@mail.fhsu.edu", logger=log)
    nucleotide_pubmed_data = pd.read_csv("DATA_Nucleotide_Pubmed_data_combined.csv")
    nucleotide_pubmed_data = nucleotide_pubmed_data.rename(columns={"European Nucleotide Archive": "ENA", "SRA Code": "SRA_Code"})
    nucleotide_data_list = []
    for row in nucleotide_pubmed_data.itertuples():
        nucleotide_item = {
            "AccessionID": row.AccessionID,
            "BioProject": row.BioProject
        }
        nucleotide_info = nucleotide_interact.fetch_nucleotide_info(row.AccessionID)
        if(len(nucleotide_info) > 0 and "GBSeq_definition" in nucleotide_info[0]):
            nucleotide_item["Title"] = nucleotide_info[0]["GBSeq_definition"]
        source_info = nucleotide_interact.extract_source_info_from_nucleotide_info(nucleotide_info)
        nucleotide_item["Source"] = source_info
        nucleotide_data_list.append(nucleotide_item)
        
    output_file_name = "DATA_nucleotide_data_list_new.json"
    with open(output_file_name, "w") as file:
        json.dump(nucleotide_data_list, file, indent=4)
    log.info(f"Data written to {output_file_name}")
    log.info("Finished.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Author|Version: '+__version__)
    # parser.add_argument("--mail", "-m", type=str, required=True,
    #                     help="Your email address (needed for querying NCBI PubMed via Entrez)")
    
    # parser.add_argument("--verbose", "-v", action="store_true", required=False, 
    #                     default=True, help="(Optional) Enable verbose logging")
    # parser.add_argument("--filepath", "-f",  type=str, required=True, 
    #                     help="(Required) Filename which contain title")
    args = parser.parse_args()
    main(args= args)