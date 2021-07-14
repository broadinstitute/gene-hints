from bs4 import BeautifulSoup
import csv
import requests
from time import perf_counter

gene_list_url = "https://raw.githubusercontent.com/eweitz/ideogram/master/data/annotations/gene-cache/homo-sapiens.tsv"
wiki_base_url = "https://en.wikipedia.org/wiki/"
custom_disambiguation = { # override for pages that don't automatically lead to to the disambiguation
    "BTD": "PageDoesNotExist", # no page for the gene, but BTD has it's own page
    "BOC": "BOC_(gene)",
    "CBS": "Cystathionine_beta_synthase",
    "DST": "Dystonin",
    "MCU": "PageDoesNotExist", # no page for the gene, but MCU has it's own page
    "MDK": "Midkine",
    "MVP": "Major_vault_protein",
    "NRL": "NRL_(gene)",
    "PDF": "PDF_(gene)",
    "RDX": "Radixin",
    "SPARC": "Osteonectin",
    "T": "Brachyury",
    "VIP": "Vasoactive_intestinal_peptide" }


# Get a list of gene symbols by downloading and processing the list of all gene symbols.
def get_gene_symbols():
    genes = []
    with requests.Session() as s:
        download = s.get(gene_list_url)
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

def get_page_name_from_symbol(symbol):
    r = requests.get(wiki_base_url + symbol)
    soup = BeautifulSoup(r.content, 'html.parser') 
    return soup.find(id="firstHeading").text.replace(" ", "_")



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
                page_name = get_page_name_from_symbol(symbol)
                # if it finds a disambiguation page, try again with the specifier
                if "_(disambiguation)" in page_name:
                    page_name = get_page_name_from_symbol(symbol + "_(gene)")
            page_names_to_symbols[page_name] = symbol
            f.write("%s\t%s\n"%(page_name, symbol))
            i+=1
            if i % 100 == 0:
                print("processed", i, "in", perf_counter() - start_time, "seconds.")

symbols = get_gene_symbols()
save_page_name_map(symbols)