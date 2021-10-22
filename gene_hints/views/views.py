"""Output TSV of Wikipedia page views for all human genes in last 24 hours

This module downloads Wikipedia view counts for a pre-computed list of pages
(English Wikipedia articles) that have been mapped to gene symbols.
For example, the page https://en.wikipedia.org/wiki/Tumor_necrosis_factor maps
to the gene symbol TNF.  The map is made by `generate_gene_page_map.py`.
"""

import argparse
import csv
from datetime import datetime, timedelta
import gzip
import os
import requests
import sys
from time import perf_counter

# TODO: Consider Dask, pandas, or NumPy to speed up TSV processing
# https://medium.com/featurepreneur/pandas-vs-dask-the-power-of-parallel-computing-994a202a74bd

# Enable importing local modules when directly calling as script
if __name__ == "__main__":
    cur_dir = os.path.join(os.path.dirname(__file__))
    sys.path.append(cur_dir + "/..")

from lib import download_gzip

class Views:

    def __init__(
            self,
            cache=0,
            output_dir="./data/",
        ):
        """Define relevant URLs and directories, do other setup
        """
        downloads_dir = output_dir + "tmp/views/"
        self.name_map_tsv_path =  downloads_dir + "gene_page_map.tsv"
        self.pageviews_path = downloads_dir + "pageviews_{time}_{count}.gz"
        self.cache = cache

        # Ensure needed directory exist
        if not os.path.exists(downloads_dir):
            os.makedirs(downloads_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        self.output_dir = output_dir
        self.output_path =\
            self.output_dir + "homo-sapiens-wikipedia-views.tsv"
        self.prev_output_path =\
            self.output_dir + "homo-sapiens-wikipedia-views-prev.tsv"

        # Override CSV field limits, to handle errors in wiki pageviews files
        # (which breaks module otherwise)
        csv.field_size_limit(sys.maxsize)

    def init_views_by_gene(self, genes_by_page):
        """Initialize map of page views by gene symbol
        """
        views_by_gene = {}
        for gene in genes_by_page.values():
            views_by_gene[gene] = 0
        print("\tFound", len(views_by_gene), "genes.")
        return views_by_gene

    def load_page_to_gene_map(self):
        """Load a map from Wikipedia page names to gene symbols from file.
        Map is made by `generate_gene_page_map.py`.
        """
        name_map = {}
        path = self.name_map_tsv_path
        print(f"Loading name map from {path}")
        with open(path, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            line_count = 0
            print("Processing file contents...")
            for row in reader:
                if line_count > 0:
                    name_map[row[0].lower()]= row[1]
                line_count += 1
        return name_map

    def get_pageviews_download_url(self, time):
        directory = time.strftime("/%Y/%Y-%m/") # format like /2021/2021-07/
        # format like pageviews-20210712-130000.gz
        filename = time.strftime("pageviews-%Y%m%d-%H0000.gz")
        base_url = "https://dumps.wikimedia.org/other/pageviews"
        return base_url + directory + filename

    def get_times_and_path(self, day, hour):
        times = {}
        # Compute hourly views filename and URL
        hours = hour + (day * 24)

        # Get last hour of yesterday in UTC time.
        # This aligns with how Wikipedia reports its pageviews, and enables
        # easily verifying counts we show to users.
        #
        # For example, if you want to verify that the counts you show for
        # human PTEN are correct, compare time shown in Gene Hints UI with:
        # https://pageviews.toolforge.org/?pages=PTEN_(gene)
        today = datetime.utcnow().date()
        yday = today - timedelta(days=1)
        yesterday_first_hour = datetime(yday.year, yday.month, yday.day, 0)

        # Using 23:00 yesterday as the baseline, iterate back in time hourly
        views_datetime = yesterday_first_hour + timedelta(hours=hours)

        times = {
            # E.g. 20211014-130000
            "machine": views_datetime.strftime("%Y%m%d-%H0000"),

            # E.g. 10/14/2021 13:00
            "human": views_datetime.strftime("%m/%d/%Y, %H:00")
        }

        path = self.pageviews_path.format(time=times["machine"], count=hours)

        return times, path, views_datetime

    def download_views_file(self, day, hour):
        """Download and save Wikipedia views dump file
        """
        times, path, views_datetime = self.get_times_and_path(day, hour)
        print("Processing the views file from", times["human"])

        # Download the file
        url = self.get_pageviews_download_url(views_datetime)
        print(f"\tDownloading Wikipedia views hourly data from {url}")

        download_gzip(url, path, cache=self.cache)

    def update_views(self, views_by_gene, row, genes_by_page):
        """Check a given row from Wikipedia pageview dump file, and add any
        views to the running total.  Takes in a dictionary of pageview counts
        per gene and a row of the pageview dump file of the format
        ["language", "page", "hourly_view_count", "always_zero_val"]
        E.g.: ["aa", "Main_Page", "4", "0"]
        """
        # log and ignore malformed rows (they happen rarely for unknonwn reasons)
        if len(row) < 4:
            print("\tEncountered malformed row:", row)
        # process the row
        else:
            page = row[1].lower()
            page_in_en_wp = row[0] in ["en", "en.m"] # Desktop or mobile
            if page_in_en_wp and page in genes_by_page:
                gene = genes_by_page[page]
                views = int(row[2])
                views_by_gene[gene] += views
        return views_by_gene

    def process_views_file(self, views_by_gene, genes_by_page, day, hour):
        """Process the downloaded and zipped views file by adding all
        relevant views to the total count.
        """
        path = self.get_times_and_path(day, hour)[1]
        print(f"\tProcessing pageview file contents at {path}")
        with open(path, "r") as f:
            reader = csv.reader(f, delimiter=" ")
            line_count = 0
            for row in reader:
                line_count += 1
                views_by_gene = self.update_views(
                    views_by_gene, row, genes_by_page
                )
                if line_count % 1000000 == 0:
                    million_lines = str(int(line_count / 1000000))
                    print("\t* Processed " + million_lines + " million lines")
        return views_by_gene

    def save_to_file(self, views_by_gene, days_ago):
        """Read the existing TSV file, and update the gene counts
        The file rows should be of the format:
        ["gene", "views", "prev_views"]
        """

        # Order the gene counts
        ordered_counts = sorted(
            views_by_gene.items(), key=lambda x: x[1], reverse=True
        )
        print("Top viewed gene pages:", dict(ordered_counts[:10]))

        # Read existing data
        prev_gene_views = {}
        prev_gene_ranks = {}

        if days_ago == 0:
            with open(self.prev_output_path, "rt") as f:
                reader = csv.reader(f, delimiter="\t")
                line_count = 0
                for row in reader:
                    if line_count > 0:
                        gene = row[0]
                        views = int(row[1])
                        prev_gene_views[gene] = views
                        prev_gene_ranks[gene] = line_count
                    line_count += 1

        if days_ago == 1:
            path = self.prev_output_path
        else:
            path = self.output_path

        # Overwrite the file with new data
        with open(path, "w") as f:
            if days_ago == 1:
                headers = [
                    "# gene",
                    "views",
                    "view_rank"
                ]
            else:
                headers = [
                    "# gene",
                    "views",
                    "view_delta",
                    "view_rank",
                    "view_rank_delta"
                ]
            f.write("\t".join(headers) + "\n")
            rank = 1
            for gene, views in ordered_counts:
                if days_ago == 1:
                    columns = [gene, views, rank]
                else:
                    view_delta = views - prev_gene_views.get(gene, 0)
                    # delta is 0 if the record did not exist before
                    rank_delta = rank - prev_gene_ranks.get(gene, rank)
                    columns = [gene, views, view_delta, rank, rank_delta]
                columns = [str(col) for col in columns]
                f.write("\t".join(columns) + "\n")
                rank += 1

        print(f"Wrote Wikipedia views output file to {path}")

    def run(self, sort_by="count", debug=False):
        """Output TSV of recent Wikipedia page views for all human genes
        """
        start_time = perf_counter()

        # Get fast but incomplete data if debugging / developing
        num_hours = 24 if not debug else 2

        genes_by_page = self.load_page_to_gene_map()

        for days_ago in range(2):
            views_by_gene = self.init_views_by_gene(genes_by_page)
            for hour in range(num_hours):
                self.download_views_file(days_ago, hour)
                views_by_gene = self.process_views_file(
                    views_by_gene, genes_by_page, days_ago, hour
                )

            self.save_to_file(views_by_gene, days_ago)

        print("Finished in", int(perf_counter() - start_time), "seconds.")

# Command-line handler
if __name__ == "__main__":

    # Output docs atop this file upon invoking --help via CLI
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    ).parse_args()
    parser.add_argument(
        "--sort-by",
        help="Metric by which to sort Wikipedia views.  Count is views.",
        choices=["count", "delta", "rank", "rank_delta"],
        default="count"
    )
    parser.add_argument(
        "--debug",
        help="Get fast but incomplete data.  Useful when developing.",
        action="store_true"
    )
    args = parser.parse_args()
    sort_by = args.sort_by
    debug = args.debug

    # Run everything!
    Views().run(sort_by, debug)
