"""Output TSV of gene citation counts over time, and associated metrics

This script queries the NCBI Entrez Gene database for citation lists for
each day in the specified timeframe.

Inspired by https://github.com/pkerpedjiev/gene-citation-counts/blob/master/README
"""

import csv
from datetime import datetime, timedelta, timezone
import gzip
import os
import requests
import subprocess
import glob

from lib import read_organisms
from enrich_citations import enrich_citations
from pmids_by_date import pmids_by_date

# Number of days in timeframe
num_days = 5

def format_date(days_before=None):
    now = datetime.now(timezone.utc)
    if (days_before):
        return (now - timedelta(days=days_before)).strftime("%Y/%m/%d")
    else:
        return now.strftime("%Y/%m/%d")

def combine_daily_pmids(output_path, daily_pmid_dir):
    """Aggregate per-day files into one file, to ease downstream processing
    """
    pmids = []

    for fp in glob.glob(daily_pmid_dir + "/*tsv"):
        with open(fp) as fd:
            rd = csv.reader(fd, delimiter="\t")
            for row in rd:
                year = row[0] # TODO: Remove this column, use filename date
                pmid = row[1] # PubMed ID, i.e. citation ID
                pmids.append(year + "\t" + pmid)

    with open(output_path, "w") as f:
        lines = "\n".join(pmids)
        f.write(lines)

def download_gzip(url, output_path):
    """Download remote gzip file, decompress, write to output path
    """
    response = requests.get(url)

    try:
        # Human-readable text, to ease debugging
        content = gzip.decompress(response.content).decode()
    except gzip.BadGzipFile as e:
        print('URL did not respond with a gzipped file: ' + url)
        raise(e)

    with open(output_path, "w") as f:
        f.write(content)

def split_ncbi_file_by_org(input_path, output_filename):
    """Split a multi-organism file from NCBI into organism-specific files

    Input file must be a TSV file with taxid as first column
    """
    org_names_by_taxid = {}
    for org in organisms:
        taxid = org["taxid"]
        name = org["scientific_name"]
        org_names_by_taxid[taxid] = name

    with open(input_path, 'r') as f:
        lines_by_org = {}
        for line in f:
            line = line.strip()
            taxid = line.split('\t')[0]
            if taxid in org_names_by_taxid:
                org = org_names_by_taxid[taxid] # scientific name
                if org in lines_by_org:
                    lines_by_org[org].append(line)
                else:
                    lines_by_org[org] = [line]

        for org in lines_by_org:
            lines = lines_by_org[org]
            with open(data_dir + org + "/" + output_filename, "w") as f:
                f.write("\n".join(lines))

cites_dir = './pubmed_citations/'
data_dir = cites_dir + 'data/'
tmp_dir = data_dir + 'tmp/'

if __name__ == "__main__":

    # Make tmp_dir, and thereby also the other dirs
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    now = datetime.now(timezone.utc)

    num_days_doubled = num_days * 2

    start_date = format_date(num_days) # E.g. 60 days ago
    end_date = format_date() # Today
    prev_start_date = format_date(num_days * 2) # E.g. 120 days ago
    prev_end_date = start_date # E.g. 60 days ago

    output_dir = tmp_dir + "timeframe"
    prev_output_dir= tmp_dir + "prev_timeframe"

    pmids_by_date(start_date, end_date, output_dir)
    pmids_by_date(prev_start_date, prev_end_date, output_dir)

    # Consolidate the publications from each day into one complete list
    pmid_dates_path = data_dir + "pmid_dates.tsv"
    prev_pmid_dates_path = data_dir + "prev_pmid_dates.tsv"

    print("Combine daily publication counts")
    combine_daily_pmids(pmid_dates_path, output_dir)
    combine_daily_pmids(prev_pmid_dates_path, prev_output_dir)

    organisms = read_organisms()

    for org in organisms:
        dir = data_dir + org["scientific_name"]
        if not os.path.exists(dir):
            os.makedirs(dir)

    print("Download gene2pubmed")
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
    output_filename = "gene2pubmed"
    output_path = data_dir + output_filename
    download_gzip(url, output_path)

    print("Split gene2pubmed by organism")
    split_ncbi_file_by_org(output_path, output_filename)

    print("Download UCSC genome annotations for each organism")
    for org in organisms:
        org_name = org["scientific_name"]
        taxid = org["taxid"]
        asm_ucsc_name = org["genome_assembly_ucsc_name"]

        # Construct parameters for fetching UCSC genome annotation files,
        # e.g. https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
        base_url = "https://hgdownload.soe.ucsc.edu/goldenPath/"
        leaf = "refGene.gtf"
        if org_name in ["homo-sapiens", "rattus-norvegicus", "felis-catus"]:
            leaf = asm_ucsc_name + "." + leaf
        url = base_url + asm_ucsc_name + "/bigZips/genes/" + leaf + ".gz"
        output_path = data_dir + org_name + "/" + leaf

        download_gzip(url, output_path)

    print("Download gene_info")
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
    output_filename = "gene_info"
    output_path = data_dir + output_filename
    download_gzip(url, output_path)

    print("Split gene_info by organism")
    split_ncbi_file_by_org(output_path, output_filename)

    # Combine genomic and citation data, compute statistics, and write to TSV
    enrich_citations(pmid_dates_path, prev_pmid_dates_path, num_days)


