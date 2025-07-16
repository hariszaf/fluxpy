import os, re
import json
import shutil
import subprocess
import requests
from datetime import datetime
from pathlib import Path

def download_compounds_and_reactions(outdir):
    """
    Downloads all .tsv files from the ModelSEED/ModelSEEDDatabase repository's
    Biochemistry directory that start with 'compounds'.
    
    outdir: Path
    """
    # GitHub repository details
    owner  = "ModelSEED"
    repo   = "ModelSEEDDatabase"
    branch = "dev"
    path   = "Biochemistry"

    # GitHub API URL to list contents of the Biochemistry directory
    api_url = f"https://api.github.com/repos/{owner}/{repo}/contents/{path}?ref={branch}"

    try:
        # Send GET request to GitHub API
        response = requests.get(api_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
    except requests.exceptions.RequestException as e:
        print(f"Error fetching repository contents: {e}")
        return

    # Parse JSON response
    contents = response.json()

    # Filter files that start with 'compound' and end with '.tsv'
    compound_files = [
        file for file in contents
        if file['type'] == 'file' and file['name'].startswith('compound') and file['name'].endswith('.tsv')
    ]

    # Create a directory to save the downloaded files
    compound_dir = outdir / "downloaded_compounds"
    os.makedirs(compound_dir, exist_ok=True)
    print("All compound files have been downloaded successfully.")

    # Download files
    download(target_files=compound_files, output_dir=compound_dir)
    print("\n\n==>Downloading compound files compoleted\n\n")

    # Filter files that start with 'reaction' and end with '.tsv'
    reaction_files = [
        file for file in contents
        if file['type'] == 'file' and file['name'].startswith('reaction') and file['name'].endswith('.tsv')
    ]

    # Create a directory to save the downloaded files
    reaction_dir = outdir / "downloaded_reactions"
    os.makedirs(reaction_dir, exist_ok=True)

    # Download files
    download(target_files=reaction_files, output_dir=reaction_dir)
    print("\n\n==> Downloading reaction files completed.")

    # Merge em
    merge_files(outdir=outdir)

    # Build mapping file
    mseed_maps_to_kegg(outdir=outdir)

    # Remove download files
    shutil.rmtree(compound_dir)
    shutil.rmtree(reaction_dir)


def download(target_files, output_dir):

    # Download each target file
    for file in target_files:
        download_url = file['download_url']
        file_name    = file['name']
        file_path    = os.path.join(output_dir, file_name)
        if os.path.exists(file_path):
            print("File already retrieved.")
            continue
        print(f"Downloading {file_name} out of {len(target_files)}...")
        try:
            file_response = requests.get(download_url)
            file_response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {file_name}: {e}")
            continue

        # Write the content to a file
        with open(file_path, 'wb') as f:
            f.write(file_response.content)

def merge_files(path_to_mseed_compounds=None, path_to_mseed_reactions=None, outdir=None):
    """
    Merges the partial MSEED files to a single one for compounds and reactions, removing the first line (header) and keeping it only once.
    """
    if path_to_mseed_compounds is None:
        path_to_mseed_compounds = outdir / "downloaded_compounds"
    if path_to_mseed_reactions is None:
        path_to_mseed_reactions = outdir / "downloaded_reactions"

    if os.path.isdir(path_to_mseed_compounds):
        print("\n==>Merging compounds..")
        subprocess.run(
            f'( head -n 1 {path_to_mseed_compounds}/compound_00.tsv ; '
            f'for file in {path_to_mseed_compounds}/*; do '
            f'    tail -n +2 "$file" ; '
            f'done ) > {outdir}/MSEED_COMPOUNDS.tsv',
            shell=True,
            check=True
        )
    if os.path.isdir(path_to_mseed_reactions):
        print("\n==>Merging reactions..")
        subprocess.run(
            f'( head -1 {path_to_mseed_reactions}/reaction_00.tsv; '
            f'for file in {path_to_mseed_reactions}/*; do tail -n +2 "$file"; done ) > {outdir}/MSEED_REACTIONS.tsv',
            shell=True,
            check=True
        )


def mseed_maps_to_kegg(mseed_compound_file=None, mseed_reaction_file=None, outdir=None):
    """
    Use the complete info of the MSEED files to come up with a 3-col tsv file
    linking ModelSEED compounds to BiGG and KEGG ones.
    """
    if mseed_compound_file is None:
        mseed_compound_file = outdir / "MSEED_COMPOUNDS.tsv"
    if mseed_reaction_file is None:
        mseed_reaction_file = outdir / "MSEED_REACTIONS.tsv"

    with open(mseed_compound_file) as f:
        mseed_compounds =  f.readlines()

    mseed_maps = {}
    # Skip header
    for line in mseed_compounds[1:]:
        parts = line.split("\t")
        mid, name, formula, mass, charge, alias = parts[0], parts[2], parts[3], parts[4], parts[7], parts[18]
        bigg = extract_pattern(alias, "BiGG")  # aliases.split("BiGG:")[1].split(";")[0].strip()
        kegg = extract_pattern(alias, "KEGG")  # aliases.split("KEGG:")[1].split(";")[0].strip()
        metanetx = extract_pattern(alias, "metanetx.chemical")  # aliases.split("KEGG:")[1].split(";")[0].strip()

        mseed_maps[mid] = {}
        if bigg is not None:
            mseed_maps[mid]['bigg'] = []
            for met in bigg:
                mseed_maps[mid]['bigg'].append(met)
        if kegg is not None:
            mseed_maps[mid]['kegg'] = []
            for ko in kegg:
                mseed_maps[mid]['kegg'].append(ko)
        if metanetx is not None:
            mseed_maps[mid]['metanetx'] = []
            for term in metanetx:
                mseed_maps[mid]['metanetx'].append(term)

        mseed_maps[mid]['name']    = name
        mseed_maps[mid]['alias']   = [x.strip() for x in alias.split("|")[0].split("Name:")[-1].split(";")]
        mseed_maps[mid]['formula'] = formula
        mseed_maps[mid]['mass']    = mass
        mseed_maps[mid]['charge']  = charge

    with open(outdir / "MSEED_COMPOUNDS.json", "w") as f:
        json.dump(mseed_maps, f)

    return mseed_maps


def extract_pattern(text, pattern):
    """
    Extracts the substring after a given pattern followed by ':' until the first '|' or tab.

    :param text: The input string.
    :param pattern: The pattern to match (e.g., 'AraCyc').
    :return: The extracted substring or None if no match is found.
    """
    regex = rf"{re.escape(pattern)}:\s*([^|\t]+)"
    match = re.search(regex, text)
    if match:
        return match.group(1).split(";")
    return None


if __name__ == "__main__":
    root  = Path(os.path.realpath(__file__)).parent
    download_compounds_and_reactions(outdir=root)
    print("\nModelSEED based files were retrieved and fixed in a fluxpy-friendly format.")

