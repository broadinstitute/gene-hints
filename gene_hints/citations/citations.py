"""Output TSV of gene citation counts over time, and associated metrics

This module downloads citation and genomic data from PubMed and other
bioinformatics databases, combines that data, computes statistics, and writes
it to organism-specific TSV files.  The output gene citation hints files
are combined with gene view hints files.  These can then help power various
applications, like the gene hints ideogram at https://broad.io/gene-hints.
"""

# Inspired by https://github.com/pkerpedjiev/gene-citation-counts/blob/master/README

import argparse
import csv
from datetime import datetime, timedelta, timezone
import gzip
import os
import requests
import glob
import sys

if __name__ == "__main__":
    sys.path.append('../..')

from lib import read_organisms
from enrich_citations import EnrichCitations
from pmids_by_date import pmids_by_date

def format_date(days_before=None):
    """Get date strings in YYYY/MM/DD format, as expected by NCBI E-utils
    """
    now = datetime.now(timezone.utc)
    if (days_before):
        return (now - timedelta(days=days_before)).strftime("%Y/%m/%d")
    else:
        return now.strftime("%Y/%m/%d")

def download_gzip(url, output_path):
    """Download remote gzip file, decompress, write to output path
    """
    response = requests.get(url)

    try:
        # Human-readable text, to ease debugging
        content = gzip.decompress(response.content).decode()
    except gzip.BadGzipFile as e:
        print("URL did not respond with a gzipped file: " + url)
        raise(e)

    with open(output_path, "w") as f:
        f.write(content)

def merge_daily_pmids(output_path, daily_pmid_dir):
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

class Citations():

    def __init__(
        self,
        cites_dir="./pubmed_citations/"
    ):
        self.cites_dir = cites_dir
        self.data_dir = cites_dir + "data/"
        self.tmp_dir = self.data_dir + "tmp/"

    def split_ncbi_file_by_org(self, input_path, output_filename, organisms):
        """Split a multi-organism file from NCBI into organism-specific files

        Input file must be a TSV file with taxid as first column
        """
        org_names_by_taxid = {}
        for org in organisms:
            taxid = org["taxid"]
            name = org["scientific_name"]
            org_names_by_taxid[taxid] = name

        with open(input_path, "r") as f:
            lines_by_org = {}
            for line in f:
                line = line.strip()
                taxid = line.split("\t")[0]
                if taxid in org_names_by_taxid:
                    org = org_names_by_taxid[taxid] # scientific name
                    if org in lines_by_org:
                        lines_by_org[org].append(line)
                    else:
                        lines_by_org[org] = [line]

            for org in lines_by_org:
                lines = lines_by_org[org]
                output_path = self.data_dir + org + "/" + output_filename
                with open(output_path, "w") as f:
                    f.write("\n".join(lines))

    def fetch_all_publications_over_time(self, path, prev_path, num_days):
        """Download IDs for articles published in the last `num_days`
        """
        start_date = format_date(num_days) # E.g. 60 days ago
        end_date = format_date() # Today
        prev_start_date = format_date(num_days * 2) # E.g. 120 days ago
        prev_end_date = start_date # E.g. 60 days ago

        output_dir = self.tmp_dir + "timeframe"
        prev_output_dir= self.tmp_dir + "prev_timeframe"

        pmids_by_date(start_date, end_date, output_dir)
        pmids_by_date(prev_start_date, prev_end_date, prev_output_dir)

        print("Combine daily publication counts")
        merge_daily_pmids(path, output_dir)
        merge_daily_pmids(prev_path, prev_output_dir)

    def fetch_all_publications_per_organism(self, organisms):
        """Get IDs for articles published about our organisms of interest
        """
        for org in organisms:
            dir = self.data_dir + org["scientific_name"]
            if not os.path.exists(dir):
                os.makedirs(dir)

        print("Download gene2pubmed")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
        output_filename = "gene2pubmed"
        output_path = self.data_dir + output_filename
        download_gzip(url, output_path)

        print("Split gene2pubmed by organism")
        self.split_ncbi_file_by_org(output_path, output_filename, organisms)

    def fetch_genome_annotations(self, organisms):
        """Download GTF files from UCSC; these contain genomic coordinates, etc.
        """
        print("Download UCSC genome annotations for each organism")
        for org in organisms:
            org_name = org["scientific_name"]
            asm_ucsc_name = org["genome_assembly_ucsc_name"]

            # Construct parameters for fetching UCSC genome annotation files,
            # e.g. https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
            base_url = "https://hgdownload.soe.ucsc.edu/goldenPath/"
            leaf = "refGene.gtf"
            if org_name in ["homo-sapiens", "rattus-norvegicus", "felis-catus"]:
                leaf = asm_ucsc_name + "." + leaf
            url = base_url + asm_ucsc_name + "/bigZips/genes/" + leaf + ".gz"
            output_path = self.data_dir + org_name + "/" + leaf

            download_gzip(url, output_path)

    def fetch_gene_info(self, organisms):
        print("Download gene_info")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
        output_filename = "gene_info"
        output_path = self.data_dir + output_filename
        download_gzip(url, output_path)

        print("Split gene_info by organism")
        self.split_ncbi_file_by_org(output_path, output_filename, organisms)

    def download_data(self, pmid_dates_path, prev_pmid_dates_path, num_days):
        """Download citation and genomic data, preparing it for enrichment
        """

        # Make tmp_dir, and thereby also the other dirs
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        # Download IDs for articles published in the last `num_days`
        self.fetch_all_publications_over_time(
            pmid_dates_path, prev_pmid_dates_path, num_days
        )

        organisms = read_organisms()

        # Download IDs for articles published about our organisms of interest.
        # We filter and join the "publications over time" and "publications per
        # organism" lists in `enrich_citations`.
        self.fetch_all_publications_per_organism(organisms)

        # Download GTF files from UCSC; these contain genomic coordinates, etc.
        self.fetch_genome_annotations(organisms)

        # Download more genomic information
        # TODO: Is data parsed from gene_info available in UCSC GTF files?
        self.fetch_gene_info(organisms)

    def run(self, num_days):
        """Output TSV of gene citation counts and related metrics over `num_days`
        """

        pmid_dates_path = self.data_dir + "pmid_dates.tsv"
        prev_pmid_dates_path = self.data_dir + "prev_pmid_dates.tsv"

        # Download citation and genomic data
        self.download_data(pmid_dates_path, prev_pmid_dates_path, num_days)

        # Combine that downloaded data, compute statistics, and write to TSV
        EnrichCitations().run(pmid_dates_path, prev_pmid_dates_path, num_days)

# Command-line handler
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--num-days",
        type=int,
        help="Number of days to analyze",
        default=180
    )
    args = parser.parse_args()
    num_days = args.num_days

    Citations().run(num_days)
