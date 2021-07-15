from bs4 import BeautifulSoup
import csv
import requests
from time import perf_counter

gene_symbol_master_list_url = "https://raw.githubusercontent.com/eweitz/ideogram/master/data/annotations/gene-cache/homo-sapiens.tsv"
wiki_base_url = "https://en.wikipedia.org/wiki/"
custom_disambiguation = { # override for pages that don't automatically lead to to the disambiguation
    "BTD": "PageDoesNotExist", # no page for the gene, but BTD has it's own page (should be fixed now!)
    "BOC": "BOC_(gene)", # Need to handle this type of disambiguation automatically in the future
    "CBS": "Cystathionine_beta_synthase",
    "CTSH": "Cathepsin_H",
    "DST": "Dystonin",
    "ESPN": "Espin_(protein)",
    "GBA": "Glucocerebrosidase",
    "GIF": "Intrinsic_factor",
    "MCU": "PageDoesNotExist", # no page for the gene, but MCU has it's own page
    "MDK": "Midkine",
    "MVP": "Major_vault_protein",
    "NRL": "NRL_(gene)",
    "PDF": "PDF_(gene)",
    "REST": "RE1-silencing_transcription_factor",
    "RFK": "PageDoesNotExist", # no page for the gene, but RFK has it's own page
    "RDX": "Radixin",
    "SPARC": "Osteonectin",
    "T": "Brachyury",
    "VIP": "Vasoactive_intestinal_peptide" }


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
def save_page_name_map(gene_symbols):
    page_names_to_symbols = {}
    i = 0
    start_time = perf_counter()
    with open('./wikipedia_trends/gene_page_map.tsv', 'w') as f:
        f.write("tpage_name\tgene_symbol\n")
        for symbol in gene_symbols:
            # handle any pages like "PDF" where it leads to an incorrect but non-ambiguous page
            if custom_disambiguation.get(symbol) is not None:
                page_name = custom_disambiguation[symbol]
            else: 
                # load the page, and see what the title is 
                page_name = get_correct_page_name(symbol)
            page_names_to_symbols[page_name] = symbol
            f.write("%s\t%s\n"%(page_name, symbol))
            i+=1
            if i % 100 == 0:
                print("processed", i, "in", perf_counter() - start_time, "seconds.")

symbols = get_gene_symbols()
save_page_name_map(symbols)