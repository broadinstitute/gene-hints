# Create gene citation hints

## Background
Gene Hints uses citations and view data for genes to enhance discoverability for gene search.  This document discusses how gene citation hint data is generated.

The citation code uses the PubMed database to collect all publications over a given number of days. Next, NCBI data is used to match the genes to publication that mention ("cite") them.  UCSC and NCBI data are then used to match gene IDs to genomic information like each gene's chromosome and coordinates. Lastly, related citation metrics are calculated and gene citation hint TSV files are output for each organism.

This approach was pioneered in [gene citation counts by Peter Kerpedjiev](https://github.com/pkerpedjiev/gene-citation-counts), and later reimplemented and extended by a team of engineers at the Broad Institute.

After these cite hint TSVs are generated, they are combined with view hint TSVs.  The combined data is then fetched by the gene hints ideogram for display and exploration in genomic context.

## Install
Create a Python 3 virtual environment, and install dependencies:
```
cd gene-hints
python3 -m venv env --copies
source env/bin/activate
pip install -r requirements.txt
```

## Run
Run `citations.py` from the repo's root directory.

```
python3 gene_hints/citations/citations.py --num-days 180
```

The command above covers almost a year: from today until 180 days ago, and 180 days before that.  While human or mouse have gene citations in almost any timeframe, analyzing longer time periods helps pick up gene citations for less intensely studies organisms like dog or cat.  It also picks up more genes in highly-studied organisms.

## Add a new organism
Currently, this code produces TSVs for five organisms (a.k.a. taxa or species): human, mouse, rat, dog, and cat. To add another organism, add a row to `organisms.tsv` in the root directory.

###  Set taxid
Each organism has an NCBI Taxonomy ID, commonly known as a `taxid`.  The taxid for human is 9606, for example.  To get an organism's taxid, search its scientific name (e.g. "Homo sapiens") in https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi.  Add this to the `taxid` column in `organisms.tsv`.

### Set genome assembly UCSC name
Many organisms have a _genome annotation_ -- a file that indicates the organism's genes and where they are located on its genome.  Using https://hgdownload.soe.ucsc.edu/downloads.html, find your organism's `genome assembly ucsc name` ending in ".gtf.gz". After you find that name, e.g. `hg38` or `rn6`, add it to `organisms.tsv`.
