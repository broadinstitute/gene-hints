from bs4 import BeautifulSoup
import csv
import os
import requests
from time import perf_counter

gene_symbol_master_list_url = "https://raw.githubusercontent.com/eweitz/ideogram/master/data/annotations/gene-cache/homo-sapiens.tsv"
wiki_base_url = "https://en.wikipedia.org/wiki/"
output_file_location = "./wikipedia_trends/gene_page_map.tsv"
custom_disambiguation = { # override for pages that don't automatically lead to to the disambiguation
    "BTD": "PageDoesNotExist", # no page related to this gene exists, but the symbol does match an actual unrelated page name 
    "MCU": "PageDoesNotExist",
    "RFK": "PageDoesNotExist", 
    "EPYC": "PageDoesNotExist",
    "DSEL": "PageDoesNotExist",
    "CHAD": "PageDoesNotExist",
    "MT4": "PageDoesNotExist",
    "CTSH": "Cathepsin_H", # no standard disambiguation for this
    "REST": "RE1-silencing_transcription_factor",
    "MGMT": "O-6-methylguanine-DNA_methyltransferase",
    "MAVS": "Mitochondrial_antiviral-signaling_protein" }

# O-6-methylguanine-DNA_methyltransferase
append_only = True # Control whether the script overwrites existing values in the file


# Get a list of gene symbols by downloading and processing the list of all gene symbols.
def get_gene_symbols():
    genes = []
    with requests.Session() as s:
        download = s.get(gene_symbol_master_list_url)
        decoded_content = download.content.decode('utf-8')
        reader = csv.reader(decoded_content.splitlines(), delimiter='\t')
        line_count = 0
        gene_info_list = list(reader)
        for row in gene_info_list:
            line_count += 1
            if line_count > 1:
                genes.append(row[2])
    print("Found", len(genes), "genes.")
    return genes


# If the global append_only variable is false, then this should return an empty dict
# Note: this returns a gene -> pagename mapping (opposite of what we generate later), so that we can look up by gene 
def get_past_genes_mapped():
    past_mapping = {}
    if not append_only:
        return past_mapping
    try:
        print("Initializing the page map with past mapping values")
        with open(output_file_location, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            line_count = 0
            for row in reader:
                if line_count > 0:
                    past_mapping[row[1]]= row[0]
                line_count += 1
    except FileNotFoundError:
        print("Unable to find existing output file. Generating a new one from scratch.")
    print("\tUsing", len(past_mapping), "previously generated values.")
    return past_mapping


# Determine whether the page is a disambiguation, or links to a disambiguation page.
def has_disambiguation(wiki_page_soup, attempted_page_name):
    has_disambiguation_link = wiki_page_soup.find("a", {"class": "mw-disambig"}) is not None
    has_link_to_gene_page = wiki_page_soup.find("a", {"href": "/wiki/" + attempted_page_name + "_(gene)"}) is not None
    has_disambiguation_in_header = "(disambiguation)" in wiki_page_soup.find(id="firstHeading")
    has_disambiguation_box = wiki_page_soup.find("div", {"id": "disambigbox"}) is not None

    return has_disambiguation_link or has_link_to_gene_page or has_disambiguation_in_header or has_disambiguation_box


# Load the wikipedia page, and see if it redirects to a page with a different name, or has a disambiguation. 
# The purpose of this method is to avoid getting the pageviews for unrelated pages, when we intend for the gene page. 
# Ex. "/wiki/GBA" redirects to "/wiki/Game_Boy_Advance". We want to detect the disambiguation so we know 
# to go to "/wiki//wiki/GBA_(gene)" instead.
def get_correct_page_name(attempted_page_name):
    gene_suffix = "_(gene)"
    r = requests.get(wiki_base_url + attempted_page_name)
    soup = BeautifulSoup(r.content, 'html.parser')
    page_name = soup.find(id="firstHeading").text.replace(" ", "_") 

    # when there is a disambiguation, we need to look for the specific gene page
    if has_disambiguation(soup, attempted_page_name) and not (gene_suffix in attempted_page_name):
        # recursively call this with the _gene suffix
        return get_correct_page_name(attempted_page_name + gene_suffix)
    else: 
        return page_name


# Load the wikipedia page for each gene and see what the actual page name is. 
# Most gene wiki pages will load when given the gene symbol, but have a different actual page name. 
def save_page_name_map(gene_symbols, past_mapping):
    page_names_to_symbols = {}
    i = 0
    start_time = perf_counter()
    file_exists = os.path.exists(output_file_location)
    # Write to the file 
    with open(output_file_location, 'a' if append_only else 'w') as f:
        if not append_only or not file_exists:
            f.write("tpage_name\tgene_symbol\n")
        for symbol in gene_symbols:
            if past_mapping.get(symbol) is None:
                # Get the page value for this symbol
                # handle any pages like "PDF" where it leads to an incorrect but non-ambiguous page
                if custom_disambiguation.get(symbol) is not None:
                    page_name = custom_disambiguation[symbol]
                else: 
                    # load the page, and see what the title is 
                    page_name = get_correct_page_name(symbol)
                page_names_to_symbols[page_name.lower()] = symbol
                f.write("%s\t%s\n"%(page_name, symbol))
            i+=1
            if i % 100 == 0:
                print("processed", i, "in", perf_counter() - start_time, "seconds.")

symbols = get_gene_symbols()
past_mapping = get_past_genes_mapped()
save_page_name_map(symbols, past_mapping)