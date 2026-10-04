__version__ = "b_thapamagar@mail.fhsu.edu|2025-02-18"

import argparse

import pandas as pd
from Bio import Entrez, SeqIO


# Function to fetch a batch of records
def fetch_batch(start, end, nucleotide_esearch_output):
    handle = Entrez.efetch(
        db="nucleotide",
        rettype="gb",
        retmode="text",
        id=",".join(nucleotide_esearch_output["IdList"][start:end]),
    )
    batch_records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    return batch_records


def main(args):
    API_KEY = (
        args.api_key
    )  # Please extract api key from NCBI website (sign in and get access to API KEY)
    Entrez.email = (
        args.mail
    )  # Mention the email address same as which is used to sign in NCBI.
    Entrez.api_key = API_KEY
    search_term = """Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] 
                NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"""

    handle = Entrez.esearch(
        db="Nucleotide",
        term=search_term,
        usehistory="y",
        retmax=100000,  # Total count was: 62173 so, to access all UIDs, set retmax to 100000
    )

    nucleotide_esearch_output = Entrez.read(handle)

    batch_size = 1000

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
            batch_records = fetch_batch(start, end)

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

            data_records.to_csv(
                "Nucleotide_Metadata.csv", mode="a", header=True, index=False
            )
            print("Records saved successfully")
            print(
                f"Fetched batch {start // batch_size + 1}: {len(batch_records)} records"
            )

        except Exception as e:
            data_records.to_csv("Nucleotide_Metadata_on_Exception.csv", index=False)
            print("Records saved because of exception encountered")
            print(f"An error occurred: {e}")


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
