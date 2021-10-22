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
    sys.path.append("../..")

from lib import read_organisms, is_cached
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

def download_gzip(url, output_path, cache=0):
    """Download remote gzip file, decompress, write to output path
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

class Citations():

    def __init__(
        self,
        cache=0,
        cites_dir="./pubmed_citations/"
    ):
        self.cites_dir = cites_dir
        self.data_dir = cites_dir + "data/"
        self.tmp_dir = self.data_dir + "tmp/"
        self.cache = cache

    def split_ncbi_file_by_org(self, input_path, output_filename, organisms):
        """Split a multi-organism file from NCBI into organism-specific files

        Input file must be a TSV file with taxid as first column
        """

        all_cached = True

        output_paths_by_org = {}
        org_names_by_taxid = {}
        for org in organisms:
            taxid = org["taxid"]
            org = org["scientific_name"]
            org_names_by_taxid[taxid] = org

            output_path = self.data_dir + org + "/" + output_filename
            output_paths_by_org[org] = output_path

            if not is_cached(output_path, self.cache, 2):
                all_cached = False

        if all_cached:
            print(
                f"All NCBI organism files for {input_path} are cached, " +
                "so not computing any."
            )
            return

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

    def merge_daily_pmids(self, output_path, daily_pmid_dir, cache=0):
        """Aggregate per-day files into one file, to ease downstream processing
        """
        pmids = []

        if is_cached(output_path, cache, 2):
            return

        for fp in glob.glob(daily_pmid_dir + "/*tsv"):
            with open(fp) as fd:
                reader = csv.reader(fd, delimiter="\t")
                for row in reader:
                    year = row[0] # TODO: Remove this column, use filename date
                    pmid = row[1] # PubMed ID, i.e. citation ID
                    pmids.append(year + "\t" + pmid)

        with open(output_path, "w") as f:
            lines = "\n".join(pmids)
            f.write(lines)

    def fetch_all_publications_over_time(self, path, prev_path, days):
        """Download IDs for articles published in the last n `days`
        """
        start_date = format_date(days) # E.g. 60 days ago
        end_date = format_date() # Today
        prev_start_date = format_date(days * 2) # E.g. 120 days ago
        prev_end_date = start_date # E.g. 60 days ago

        output_dir = self.tmp_dir + "timeframe"
        prev_output_dir= self.tmp_dir + "prev_timeframe"

        pmids_by_date(start_date, end_date, output_dir, self.cache)
        pmids_by_date(prev_start_date, prev_end_date, prev_output_dir, self.cache)

        print("Combine daily publication counts")
        self.merge_daily_pmids(path, output_dir)
        self.merge_daily_pmids(prev_path, prev_output_dir)

    def fetch_all_publications_per_organism(self, organisms):
        """Get IDs for articles published about our organisms of interest
        """
        for org in organisms:
            dir = self.data_dir + org["scientific_name"]
            if not os.path.exists(dir):
                os.makedirs(dir)

        print("Download gene2pubmed")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
        output_name = "gene2pubmed"
        output_path = self.data_dir + output_name
        download_gzip(url, output_path, self.cache)

        print("Split gene2pubmed by organism")
        self.split_ncbi_file_by_org(output_path, output_name, organisms)

    def fetch_gene_info(self, organisms):
        print("Download gene_info")
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
        output_name = "gene_info"
        output_path = self.data_dir + output_name
        download_gzip(url, output_path, self.cache)

        print("Split gene_info by organism")
        self.split_ncbi_file_by_org(output_path, output_name, organisms)

    def download_data(self, pmid_dates_path, prev_pmid_dates_path, days):
        """Download citation and genomic data, preparing it for enrichment
        """

        # Make tmp_dir, and thereby also the other dirs
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        # Download IDs for articles published in the last `days`
        self.fetch_all_publications_over_time(
            pmid_dates_path, prev_pmid_dates_path, days
        )

        organisms = read_organisms()

        # Download IDs for articles published about our organisms of interest.
        # We filter and join the "publications over time" and "publications per
        # organism" lists in `enrich_citations`.
        self.fetch_all_publications_per_organism(organisms)

        # Download more genomic information
        # TODO: Is data parsed from gene_info available in UCSC GTF files?
        self.fetch_gene_info(organisms)

    def run(self, days, sort_by="count"):
        """Output TSV of gene citation counts and related metrics over `days`
        """

        pmid_dates_path = self.data_dir + "pmid_dates.tsv"
        prev_pmid_dates_path = self.data_dir + "prev_pmid_dates.tsv"

        # Download citation and genomic data
        self.download_data(pmid_dates_path, prev_pmid_dates_path, days)

        # Combine that downloaded data, compute statistics, and write to TSV
        EnrichCitations().run(
            pmid_dates_path, prev_pmid_dates_path, days, sort_by
        )

# Command-line handler
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--days",
        type=int,
        help="Number of days to analyze",
        default=180
    )
    parser.add_argument(
        "--sort-by",
        help="Metric by which to sort PubMed citations.  Count is citations.",
        choices=["count", "delta", "rank", "rank_delta"],
        default="count"
    )
    parser.add_argument(
        "--cache",
        help=(
            "Get fast but incomplete data.  Useful to develop.  Levels:" +
                "0: Don't cache.  " +
                "1: Cache download but not compute.  " +
                "2: like debug=1, and cache intermediate compute.  " +
                "(default: %(default)i)"
        ),
        choices=[0, 1, 2],
        default=0
    )
    args = parser.parse_args()
    days = args.days
    sort_by = args.sort_by
    cache = args.cache

    Citations(cache).run(days, sort_by)
