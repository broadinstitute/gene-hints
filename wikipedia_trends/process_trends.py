import csv
from datetime import datetime, timedelta, timezone
import gzip
import os
import requests
import sys
from time import perf_counter

# TODO consider using numpy to speed up TSV reading and writing 

# URLS and Directory locations
wiki_pageviews_base_url = "https://dumps.wikimedia.org/other/pageviews"
wiki_trends_dir = "./wikipedia_trends/"
name_map_tsv_location =  wiki_trends_dir + "gene_page_map.tsv"
downloads_dir = wiki_trends_dir + "downloads/"
pageviews_download_location = downloads_dir + "pageviews{count}.gz"
output_file_location = "./data/homo-sapiens-wikipedia-trends.tsv"


# Other global variables and settings
num_hours_to_process = 24
most_recent_datetime = datetime.now(timezone.utc) - timedelta(hours=1) # One hour before the current time in UTC

# override the csv field limits, to handle errors in wiki trends files (which will break the whole script otherwise)
csv.field_size_limit(sys.maxsize) 



# Initialize our gene counts as zero by downloading and processing the list of all gene symbols.
def init_gene_counts(page_to_gene_map):
    gene_counts = {}
    for gene in page_to_gene_map.values():
        gene_counts[gene] = 0
    print("\tFound", len(gene_counts), "genes.")
    return gene_counts


# Load a map from wikipedia page names to gene symbols from file.
# (This map is generated by running the other script).
def load_page_name_to_gene_map():
    name_map = {}
    print("Loading the name mapping from file...")
    with open(name_map_tsv_location, "rt") as f:
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
    filename = time.strftime("pageviews-%Y%m%d-%H0000.gz") # format like pageviews-20210712-130000.gz
    return wiki_pageviews_base_url + directory + filename


# Download and save the Wikipedia trends dump file 
def download_trends_file(file_num):
    # Ensure that the downloads directory exists
    if not os.path.exists(downloads_dir):
        os.makedirs(downloads_dir)
    # Generate the hourly trends filename and URL
    trends_datetime = most_recent_datetime + timedelta(hours=-file_num)
    pageviews_url = get_pageviews_download_url(trends_datetime)
    print("Processing the trends file from", trends_datetime.strftime("%m/%d/%Y, %H:00"))
    # Download the file 
    print("\tDownloading wikipedia trends hourly data...")
    with requests.Session() as s:
        response = s.get(pageviews_url)
        with open(pageviews_download_location.format(count=file_num), "wb") as f:
            f.write(response.content)


# Checks a given row from the wikipedia pageview dump file, and adds any views to the running total. 
# Takes in a dictionary of pageview counts per gene and a row of the 
# pageview dump file of the format ["language", "page_name", "hourly_view_count", "always_zero_val"]
# Ex. ["aa", "Main_Page", "4", "0"]
def add_row_to_gene_count(counts_dict, row, page_to_gene_map):
    # log and ignore malformed rows (they happen rarely for uknonwn reasons)
    if len(row) < 4:
        print("\tEncountered malformed row:", row)
    # process the row 
    else:
        page_name = row[1].lower()
        if row[0] == "en" and page_to_gene_map.get(page_name) is not None:
            gene_symbol = page_to_gene_map[page_name]
            view_count = int(row[2])
            counts_dict[gene_symbol] += view_count


# Process the downloaded and zipped trends file by adding all
# relevant views to the total count. 
def process_trends_file(gene_counts, page_to_gene_map, file_num):
    print("\tProcessing pageview file contents...")
    with gzip.open(pageviews_download_location.format(count=file_num), "rt") as f:
        reader = csv.reader(f.read().splitlines(), delimiter=" ")
        line_count = 0
        for row in reader:
            line_count += 1
            add_row_to_gene_count(gene_counts, row, page_to_gene_map)
            if line_count % 1000000 == 0:
                print("\t-- Processed", int(line_count / 1000000), " million lines. --")


# Read the existing tsv file, and update the gene counts 
# The file rows should be of the format; ["gene_symbol", "daily_page_views", "prev_daily_page_views"]
def save_counts_to_file(view_counts):
    # Order the gene counts 
    ordered_counts = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
    print("Top viewed gene pages:", dict(ordered_counts[:10]))
    # Read in the existing data
    prev_view_counts = {}
    prev_gene_ranks = {}
    print("Updating the wikipedia trends output file...")
    with open(output_file_location, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        line_count = 0
        for row in reader:
            if line_count > 0:
                gene_symbol = row[0]
                view_count = int(row[1])
                prev_view_counts[gene_symbol] = view_count
                prev_gene_ranks[gene_symbol] = line_count
            line_count += 1
    # Overwrite the file with new data
    with open(output_file_location, "w") as f:
        f.write("gene_symbol\twikipedia_daily_page_views\twikipedia_daily_page_views_change_from_previous_day\tview_rank\tview_rank_delta\n")
        rank = 1
        for gene, views, in ordered_counts:
            prev_views = prev_view_counts.get(gene, 0)
            rank_delta = prev_gene_ranks.get(gene, rank) - rank # delta is 0 if the record did not exist before
            f.write("%s\t%s\t%s\t%s\t%s\n"%(gene, views, prev_views, rank, rank_delta))
            rank += 1


# Run everything!
start_time = perf_counter()

page_to_gene_map = load_page_name_to_gene_map()
gene_counts = init_gene_counts(page_to_gene_map)

for file_num in range(num_hours_to_process):
    download_trends_file(file_num)
    process_trends_file(gene_counts, page_to_gene_map, file_num)

save_counts_to_file(gene_counts)

print("Finished in", int(perf_counter() - start_time), "seconds.")