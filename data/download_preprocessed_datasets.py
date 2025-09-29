import gdown
import os
from constants import BASE_PATH_PREPROCESSED, BASE_PATH_DATA

CANSIG_PP_URLS = {
    "raw_data": "https://drive.google.com/drive/folders/1Kp_4gvEwGMOdiEYVupyQAns_yryMMYGH?usp=sharing"
    "cansig_processed": "https://drive.google.com/drive/folders/1KITDZ1pwYDJtZSe0aJln-tn1ynpPZm4z?usp=drive_link"
}

PREPROCESSED_URLS = {
    "preprocessed": "https://drive.google.com/drive/folders/1cDg_9kVfiNynXVXxfIR8iIxPxAi8c2BN?usp=drive_link"
}

ANNOTATION_URLS = {
    "annotations": "https://drive.google.com/drive/folders/1QnGMeSJ2J1jeHF1JXuxb1yZYE6fr6GLW?usp=drive_link",
    "dgex_genes": "https://drive.google.com/drive/folders/1SYBMjmnUFCLJjciwWWAEm9zmdfUTn985?usp=drive_link",
}


def download_with_check(name, url, target_path, is_folder=False, file_extension=""):
    target_file = os.path.join(target_path, f"{name}{file_extension}")

    print(f"Downloading {name}...")

    if os.path.exists(target_file):
        print(f">> {target_file} already exists")
        return

    if is_folder:
        gdown.download_folder(url, output=target_file, quiet=False)
    else:
        gdown.download(url, target_file, quiet=False)

    print(f">> {target_file} downloaded")


if __name__ == "__main__":
    if not os.path.exists(BASE_PATH_DATA) or not os.path.isdir(BASE_PATH_PREPROCESSED):
        raise ValueError(
            f"Base path {BASE_PATH_DATA} or {BASE_PATH_PREPROCESSED} do not exist or is not a directory."
            "Please create the directories with `python constants.py`"
            )

    ## Download preprocessed datasets
    for dataset, url in PREPROCESSED_URLS.items():
        download_with_check(
            name=dataset,
            url=url,
            target_path=BASE_PATH_DATA,
            is_folder=True,
            file_extension="",
        )

    ## Download cansig processed datasets
    for annotation_folder, url in CANSIG_PP_URLS.items():
        download_with_check(
            name=annotation_folder,
            url=url,
            target_path=BASE_PATH_DATA,
            is_folder=True,
            file_extension="",
        )

    ## Download annotation folders
    for annotation_folder, url in ANNOTATION_URLS.items():
        download_with_check(
            name=annotation_folder,
            url=url,
            target_path=BASE_PATH_DATA,
            is_folder=True,
            file_extension="",
        )
