"""Output TSV of Wikipedia page views for all human genes in last 24 hours

This module downloads Wikipedia view counts for a pre-computed list of pages
(English Wikipedia articles) that have been mapped to gene symbols.
For example, the page https://en.wikipedia.org/wiki/Tumor_necrosis_factor maps
to the gene symbol TNF.  The map is made by `generate_gene_page_map.py`.
"""

import argparse
import csv
from datetime import datetime, timedelta, timezone
import gzip
import os
import requests
import sys
from time import perf_counter

# TODO consider using numpy to speed up TSV reading and writing

# URLs and directory paths
wiki_pageviews_base_url = "https://dumps.wikimedia.org/other/pageviews"
wiki_views_dir = "./gene_hints/views/"
name_map_tsv_path =  wiki_views_dir + "gene_page_map.tsv"
downloads_dir = wiki_views_dir + "downloads/"
pageviews_download_path = downloads_dir + "pageviews{count}.gz"
output_path = "./data/homo-sapiens-wikipedia-views.tsv"

# Other global variables and settings
num_hours_to_process = 24

# One hour before the current time in UTC
hour_before_present = datetime.now(timezone.utc) - timedelta(hours=2)

# override CSV field limits, to handle errors in wiki views files
# (which will break the whole script otherwise)
csv.field_size_limit(sys.maxsize)

def init_views_by_gene(genes_by_page):
    """Initialize map of page views by gene symbol
    """
    views_by_gene = {}
    for gene in genes_by_page.values():
        views_by_gene[gene] = 0
    print("\tFound", len(views_by_gene), "genes.")
    return views_by_gene

def load_page_to_gene_map():
    """Load a map from Wikipedia page names to gene symbols from file.
    (This map is generated by running the other script).
    """
    name_map = {}
    print("Loading the name mapping from file...")
    with open(name_map_tsv_path, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        line_count = 0
        print("Processing file contents...")
        for row in reader:
            if line_count > 0:
                name_map[row[0].lower()]= row[1]
            line_count += 1
    return name_map

def get_pageviews_download_url(time):
    directory = time.strftime("/%Y/%Y-%m/") # format like /2021/2021-07/
    # format like pageviews-20210712-130000.gz
    filename = time.strftime("pageviews-%Y%m%d-%H0000.gz")
    return wiki_pageviews_base_url + directory + filename

def download_views_file(hour):
    """Download and save Wikipedia views dump file
    """
    # Ensure that the downloads directory exists
    if not os.path.exists(downloads_dir):
        os.makedirs(downloads_dir)
    # Generate the hourly views filename and URL
    views_datetime = hour_before_present + timedelta(hours=-hour)
    url = get_pageviews_download_url(views_datetime)
    friendly_time = views_datetime.strftime("%m/%d/%Y, %H:00")
    print("Processing the views file from", friendly_time)
    # Download the file
    print(f"\tDownloading Wikipedia views hourly data from {url}")
    with requests.Session() as s:
        response = s.get(url)
    with open(pageviews_download_path.format(count=hour), "wb") as f:
        f.write(response.content)

def update_views(views_by_gene, row, genes_by_page):
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
        if row[0] == "en" and page in genes_by_page:
            gene = genes_by_page[page]
            views = int(row[2])
            views_by_gene[gene] += views
    return views_by_gene

def process_views_file(views_by_gene, genes_by_page, hour):
    """Process the downloaded and zipped views file by adding all
    relevant views to the total count.
    """
    path = pageviews_download_path.format(count=hour)
    print(f"\tProcessing pageview file contents at {path}")
    with gzip.open(path, "rt") as f:
        reader = csv.reader(f.read().splitlines(), delimiter=" ")
        line_count = 0
        for row in reader:
            line_count += 1
            views_by_gene = update_views(views_by_gene, row, genes_by_page)
            if line_count % 1000000 == 0:
                million_lines = str(int(line_count / 1000000))
                print("\t-- Processed " + million_lines + " million lines. --")
    return views_by_gene

def save_counts_to_file(views_by_gene):
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

    if os.path.exists(output_path):
        with open(output_path, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            line_count = 0
            for row in reader:
                if line_count > 0:
                    gene = row[0]
                    views = int(row[1])
                    prev_gene_views[gene] = views
                    prev_gene_ranks[gene] = line_count
                line_count += 1

    print("Updating Wikipedia views output file...")

    # Overwrite the file with new data
    with open(output_path, "w") as f:
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
            prev_views = prev_gene_views.get(gene, 0)
            # delta is 0 if the record did not exist before
            rank_delta = prev_gene_ranks.get(gene, rank) - rank
            columns = [gene, views, prev_views, rank, rank_delta]
            columns = [str(col) for col in columns]
            f.write("\t".join(columns) + "\n")
            rank += 1

def gene_views():
    """Output TSV of Wikipedia page views for all human genes in last 24 hours
    """
    start_time = perf_counter()

    genes_by_page = load_page_to_gene_map()
    views_by_gene = init_views_by_gene(genes_by_page)

    for hour in range(num_hours_to_process):
        download_views_file(hour)
        views_by_gene = process_views_file(views_by_gene, genes_by_page, hour)

    save_counts_to_file(views_by_gene)

    print("Finished in", int(perf_counter() - start_time), "seconds.")

# Command-line handler
if __name__ == "__main__":

    # Output docs atop this file upon invoking --help via CLI
    argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    ).parse_args()

    # Run everything!
    gene_views()
