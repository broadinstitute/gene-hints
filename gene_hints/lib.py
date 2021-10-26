"""Functions used by multiple Gene Hints modules
"""

import csv
import os
import gzip
import requests

def read_organisms():
    """Read organisms TSV and parse taxon / species data

    TODO:
    * Pull these values from NCBI EUtils API instead of hard-coding
        - Do not pull down entire NCBI Taxonomy DB.  It's much too big for
            this use case, slowing and complicating development.

    * Make `organisms` a CLI parameter
    """
    organisms = []

    with open("./organisms.tsv") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) == 0 or row[0][0] == "#":
                continue

            # Convert e.g. "Homo sapiens" to machine-friendlier "homo-sapiens"
            name = row[1]
            name = name.lower().replace(" ", "-")

            organism = {
                "common_name": row[0], # e.g. human
                "scientific_name": name, # e.g. homo-sapiens (see note above)
                "taxid": row[2] # e.g. 9606 ("NCBI Taxonomy ID")
            }

            organisms.append(organism)

    return organisms

def is_cached(path, cache, threshold):
    """Determine if file path is already available, per cache and threshold.

    `cache` level is set by pipeline user; `threshold` by the calling function.

    See `--help` CLI output for description of `cache` levels.
    """

    if cache >= threshold:
        if threshold == 1:
            action = "download"
        elif threshold == 2:
            action = "comput"

        if os.path.exists(path):
            print(f"Using cached copy of {action}ed file {path}")
            return True
        else:
            print(f"No cached copy exists, so {action}ing {path}")
            return False
    return False

def download_gzip(url, output_path, cache=0):
    """Download gzip file, decompress, write to output path; use optional cache

    Cached files can help speed development iterations by > 2x, and some
    development scenarios (e.g. on a train or otherwise without an Internet
    connection) can be impossible without it.
    """

    if is_cached(output_path, cache, 1):
        return

    response = requests.get(url)

    try:
        # Human-readable text, to ease debugging
        content = gzip.decompress(response.content).decode()
    except gzip.BadGzipFile as e:
        print("URL did not respond with a gzipped file: " + url)
        raise(e)

    with open(output_path, "w") as f:
        f.write(content)
