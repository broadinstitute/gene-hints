"""Combine TSVs of Wikipedia views and PubMed citations into one TSV file

Example usage:
python3 gene_hints/merge_hints.py --citations-path=data/homo-sapiens-pubmed-citations.tsv --views-path=data/homo-sapiens-wikipedia-views.tsv --output-path=data/homo-sapiens-gene-hints.tsv
"""

import argparse
import csv

def tsv_to_dict(tsv_path):
    print("Processing " + tsv_path)
    tsv_dict = {}
    headers = []
    with open(tsv_path) as fd:
        rd = csv.reader(fd, delimiter="\t")
        is_header_row = True
        for row in rd:
            if is_header_row: # save headers to be used later
                is_header_row = False
                headers = row
                continue
            gene = row[0]
            remaining_rows = row[1:] # everything after the gene
            tsv_dict[gene] = remaining_rows

    print("Row count: " + str(len(tsv_dict.keys())))
    first_dict_entry = tsv_dict[list(tsv_dict.keys())[0]]
    print("first_dict_entry: " + str(first_dict_entry))

    return tsv_dict, headers

def merge_dicts(tsv_dict_1, tsv_dict_2):
    print("Starting merge_dicts")

    unmatched_rows_counter = 0
    for gene in tsv_dict_1:
        if gene in tsv_dict_2.keys():
            tsv_dict_1[gene] = tsv_dict_1[gene] + tsv_dict_2[gene]
        else:
            empty_views_data = ["0", "0"] # as strings to match rest of data
            tsv_dict_1[gene] = tsv_dict_1[gene] + empty_views_data
            unmatched_rows_counter = unmatched_rows_counter + 1
    print("unmatched_rows_counter: " + str(unmatched_rows_counter))

    first_dict_entry = tsv_dict_1[list(tsv_dict_1.keys())[0]]
    print("first_dict_entry: " + str(first_dict_entry))

    return tsv_dict_1

def merge_headers(headers_1, headers_2):
    """Merge two lists of headers; assumes first column is same in both
    """
    print("Merging: " + str(headers_1))
    print("With: " + str(headers_2))
    if headers_1[0] != headers_2[0]:
        exit("Error: First column of the headers do not match. Aborting!")

    merged_headers = []
    headers_2_without_primary_key = headers_2[1:]

    merged_headers = headers_1 + headers_2_without_primary_key
    print("merged_headers list: " + str(merged_headers))
    return merged_headers

def dict_to_tsv(merged_dict, merged_headers, output_path):
    with open(output_path, "wt") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(merged_headers)
        for key in merged_dict:
            row_data = [key] + merged_dict[key]
            tsv_writer.writerow(row_data)

def merge_hints(citations_path, views_path, output_path):
    """Combine TSVs of Wikipedia views and PubMed citations into one TSV file
    """
    citations_data, citations_headers = tsv_to_dict(citations_path)
    views_data, views_headers = tsv_to_dict(views_path)

    # Citation data has more rows than views data, so use it as "base" TSV
    merged_dict = merge_dicts(citations_data, views_data)
    merged_headers = merge_headers(citations_headers, views_headers)

    dict_to_tsv(merged_dict, merged_headers, output_path)

# Command-line handler
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, # Output docs atop this file
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--citations-path", help="Path to PubMed citation hints TSV")
    parser.add_argument("--views-path", help="Path to Wikipedia view hints TSV")
    parser.add_argument("--output-path", help="Path for TSV output")
    args = parser.parse_args()

    merge_hints(args.citations_path, args.views_path, args.output_path)
