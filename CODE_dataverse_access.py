#!/usr/bin/env python3
__version__ = "b_thapamagar@mail.fhsu.edu|2025-02-18"

import argparse
import json
import requests
from pathlib import Path


def create_dataverse(server_url, api_key, dataverse_metadata_info):
    dataverse_creation_headers = (
        {"X-Dataverse-key": api_key, "Content-Type": "application/json"}
        if api_key
        else {}
    )
    resp = requests.post(
        f"{server_url}/api/dataverses/demo",
        headers=dataverse_creation_headers,
        json=dataverse_metadata_info,  # requests serialises for you
    )
    print(resp)
    return resp


# def create_dataset(server_url, api_key):
#     dataverse_creation_headers = (
#         {"X-Dataverse-key": api_key, "Content-Type": "application/json"}
#         if api_key
#         else {}
#     )
#     resp = requests.post(
#         f"{server_url}/api/dataverses/demo",
#         headers=dataverse_creation_headers,
#         json=dataverse_metadata_info,  # requests serialises for you
#     )
#     print(resp)
#     return resp


def get_dataverse_info(server_url, dataverse_alias_name, api_key):
    headers = {"X-Dataverse-key": api_key}
    dataverse_response = requests.get(
        f"{server_url}/api/dataverses/{dataverse_alias_name}", headers=headers
    )
    if dataverse_response.status_code == 200:
        print(f"✅ Dataverse '{dataverse_alias_name}' exists.")
        print("Title:", dataverse_response.json()["data"]["name"])
    elif dataverse_response.status_code == 404:
        print(f"❌ Dataverse '{dataverse_alias_name}' does NOT exist.")
    else:
        print(
            f"⚠️ Error {dataverse_response.status_code}: {dataverse_response.text[:500]}"
        )
    return dataverse_response


def get_dataset_version_info(server_url, persistent_id, api_key):
    headers = {"X-Dataverse-key": api_key}
    get_version_info_response = requests.get(
        f"{server_url}/api/datasets/:persistentId/versions",
        params={"persistentId": persistent_id},
        headers=headers,
    )
    return get_version_info_response


def publish_draft_dataset(server_url, dataset_id, api_key):
    headers = {"X-Dataverse-key": api_key}
    release_type = "major"
    publish_url = (
        f"{server_url}/api/datasets/{dataset_id}/actions/:publish?type={release_type}"
    )
    publish_response = requests.post(publish_url, headers=headers)
    return publish_response


def upload_files_to_dataset(server_url, persistent_id, api_key):
    file_path = "DataVerse/generated.json"
    headers = {"X-Dataverse-key": api_key}
    # The same JSON block you passed after -F 'jsonData=…'
    metadata = dict(
        description="Blue skies!", categories=["Lily", "Rosemary", "Jack of Hearts"]
    )
    files = {
        # form-field name → (filename, file-object)
        "file": open(file_path, "rb"),
        # second multipart part: (no filename, json string, mimetype)
        "jsonData": (None, json.dumps(metadata), "application/json"),
    }
    upload_response = requests.post(
        f"{server_url}/api/datasets/:persistentId/add",
        params={"persistentId": persistent_id},
        headers=headers,
        files=files,
    )
    print("Status:", upload_response.status_code)
    print(upload_response.json() if upload_response.ok else upload_response.text)
    return upload_response


def main(args):
    api_key = "baf3c7fd-f8aa-424a-b295-c81fd23764ca"
    server_url = "https://demo.dataverse.org"

    data_verse_metadata_filepath = Path("DataVerse/dataverse-complete.json")
    dataverse_metadata_info = json.loads(
        data_verse_metadata_filepath.read_text(encoding="utf‑8")
    )  # validates JSON syntax

    dataverse_alias_name = dataverse_metadata_info["alias"]
    dataverse_metadata_info = get_dataverse_info(
        server_url, dataverse_alias_name, api_key
    )

    headers = {"X-Dataverse-key": api_key}
    persistent_id = "doi:10.70122/FK2/JOTITY"
    # download_api = "https://demo.dataverse.org/api/access/dataset/:persistentId/?persistentId=doi:10.70122/FK2/N2XGBJ"

    # # Build & send the request
    get_dataset_response = requests.get(
        f"{server_url}/api/datasets/:persistentId/",
        params={"persistentId": persistent_id},  # requests does the URL‑encoding
        headers=headers,
    )
    print(get_dataset_response)

    dataset_id = get_dataset_response.json()["data"]["id"]

    upload_response = upload_files_to_dataset(server_url, persistent_id, api_key)
    print(upload_response)
    #### Dataverse contents (list of dataset and inner datverse)

    dv_response = requests.get(
        f"{server_url}/api/dataverses/{dataverse_alias_name}/contents",  # requests does the URL‑encoding
        headers=headers,
    )

    print(dv_response)

    ### Get Dataset Version information
    get_version_info_response = get_dataset_version_info(
        server_url, persistent_id, api_key
    )

    print(f"version Response {get_version_info_response}")

    # STEP 2 – publish the draft
    # ---------------------------------------------------------------------
    publish_response = publish_draft_dataset(server_url, dataset_id, api_key)
    print(publish_response)

    print("Completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Author|Version: " + __version__)
    # parser.add_argument("--mail", "-m", type=str, required=True,
    #                     help="Your email address (needed for querying NCBI PubMed via Entrez)")

    # parser.add_argument("--verbose", "-v", action="store_true", required=False,
    #                     default=True, help="(Optional) Enable verbose logging")
    # parser.add_argument("--filepath", "-f",  type=str, required=True,
    #                     help="(Required) Filename which contain title")
    args = parser.parse_args()
    main(args=args)
