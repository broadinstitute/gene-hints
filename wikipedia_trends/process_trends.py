import csv
from datetime import datetime
import gzip
import requests
import sys
from time import perf_counter

# TODO: iterate through downloaded files 
# TODO: abstract the local file names

gene_list_url = "https://raw.githubusercontent.com/eweitz/ideogram/master/data/annotations/gene-cache/homo-sapiens.tsv"
pageviews_base_url = "https://dumps.wikimedia.org/other/pageviews"

csv.field_size_limit(sys.maxsize) # TODO: decide whether this is needed 

def get_pageviews_download_url(time):
    directory = time.strftime("/%Y/%Y-%m/") # format like /2021/2021-07/
    filename = time.strftime("pageviews-%Y%m%d-%H0000.gz") # format like pageviews-20210712-130000.gz
    return pageviews_base_url + directory + filename

# Checks a given row from the wikipedia pageview dump file, and adds any views to the running total. 
# Takes in a dictionary of pageview counts per gene and a row of the 
# pageview dump file of the format ['language', 'page_name', 'hourly_view_count', 'always_zero_val']
# Ex. ['aa', 'Main_Page', '4', '0']
def add_row_to_gene_count(counts_dict, row):
    page_name = row[1]
    if row[0] == "en" and counts_dict.get(page_name) is not None:
        view_count = int(row[2])
        if view_count > 0:
            counts_dict[page_name] += view_count
            print("Adding", view_count, "views to", page_name)


# Initialize our gene counts as zero by downloading and processing the list of all gene symbols.
def init_gene_counts():
    gene_counts = {}
    with requests.Session() as s:
        download = s.get(gene_list_url)
        decoded_content = download.content.decode('utf-8')
        reader = csv.reader(decoded_content.splitlines(), delimiter='\t')
        line_count = 0
        gene_info_list = list(reader)
        for row in gene_info_list:
            line_count += 1
            if line_count > 1:
                gene_counts[row[2]] = 0
    print("Found", len(gene_counts), "genes.")
    return gene_counts


# Download and save the Wikipedia trends dump file 
def download_trends_file(pageviews_url):
    print("Downloading wikipedia trends hourly data from", pageviews_url)
    with requests.Session() as s:
        response = s.get(pageviews_url)
        with open('./wikipedia_trends/downloads/test.gz', 'wb') as f:
            f.write(response.content)


# Process the downloaded and zipped trends file by adding all
# relevant views to the total count. 
def process_trends_file(gene_counts):
    print("Unzipping the file...")
    with gzip.open('./wikipedia_trends/downloads/test.gz', 'rt') as f:
        reader = csv.reader(f.read().splitlines(), delimiter=' ')
        line_count = 0
        print("Processing file contents...")
        for row in reader:
            line_count += 1
            add_row_to_gene_count(gene_counts, row)
            if line_count % 1000000 == 0:
                print("-- Processed", line_count, "lines. --")


# Save to file and print the top viewed gene pages. 
def output_top_counts(gene_counts):
    # Get the top counts per gene 
    top_counts = dict(sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:10])
    print("Top viewed gene pages:", top_counts)
    with open('./wikipedia_trends/wikipedia_top_viewed_genes.tsv', 'w') as f:
        f.write("gene_symbol\tpage_views\n")
        for gene, views, in gene_counts.items():
            f.write("%s\t%s\n"%(gene, views))


# Get a list of all gene symbols for which there is a disambiguation wiki page.
# Iterate through all of the wikipedia pages in a hourly pageview file and identify 
# the ones which have a "_(disambiguation)" page.
def get_ambiguous_gene_symbols(gene_counts_dict):
    disambig_str = "_(disambiguation)" 
    ambiguous_symbols = []
    print("Unzipping the file...")
    with gzip.open('./wikipedia_trends/downloads/test.gz', 'rt') as f:
        reader = csv.reader(f.read().splitlines(), delimiter=' ')
        print("Processing file contents...")
        for row in reader:
            # check if it starts with a gene symbol and has a disambiguation
            page_name = row[1]
            if disambig_str in page_name:
                base_name = page_name.replace(disambig_str, "")
                if gene_counts_dict.get(base_name) is not None:
                    ambiguous_symbols.append(base_name)
    return ambiguous_symbols




# Let's run this thing! 
start_time = perf_counter()
gene_counts = init_gene_counts()
# pageviews_url = get_pageviews_download_url(datetime.now())
# download_trends_file(pageviews_url)
# process_trends_file(gene_counts)
# output_top_counts(gene_counts)
print(get_ambiguous_gene_symbols(gene_counts))
print("Finished in", perf_counter() - start_time, "seconds.")
