"""Functions used by multiple Gene Hints modules
"""

import csv
import os

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
        rd = csv.reader(f, delimiter="\t")
        for row in rd:
            if len(row) == 0 or row[0][0] == '#':
                continue

            # Convert e.g. "Homo sapiens" to machine-friendlier "homo-sapiens"
            name = row[1]
            name = name.lower().replace(' ', '-')

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
