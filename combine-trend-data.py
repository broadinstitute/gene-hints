""" 
Combines wikipedia view and citation count data tsvs
"""
usage = """
    python3 combine-trend-data.py data/homo-sapiens-citation-information.tsv data/homo-sapiens-wikipedia-trends.tsv data/homo-sapiens-gene-details.tsv
    """

import csv
import argparse

def tsv_to_dict(tsv):
    print("processing " + tsv)
    tsv_dict = {}
    headers = []
    with open(tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        isHeaderRow = True
        for row in rd:
            if isHeaderRow: # save headers to be used later
                isHeaderRow = False
                headers = row
                continue
            gene_symbol = row[0]
            remaining_rows = row[1:] # everything after the gene_symbol
            tsv_dict[gene_symbol] = remaining_rows

    print("row count: " + str(len(tsv_dict.keys())))
    first_dict_entry = tsv_dict[list(tsv_dict.keys())[0]]
    print("first_dict_entry: " + str(first_dict_entry))

    return tsv_dict, headers

def combine_dicts(tsv_dict_1, tsv_dict_2):
    print("Starting combine_dicts")

    unmatched_rows_counter = 0
    for gene_symbol in tsv_dict_1:
        if gene_symbol in tsv_dict_2.keys():
            tsv_dict_1[gene_symbol] = tsv_dict_1[gene_symbol] + tsv_dict_2[gene_symbol]
        else:
            empty_wiki_data = ['0', '0'] # as strings to match the rest of the data
            tsv_dict_1[gene_symbol] = tsv_dict_1[gene_symbol] + empty_wiki_data
            unmatched_rows_counter = unmatched_rows_counter + 1
    print("unmatched_rows_counter: " + str(unmatched_rows_counter))

    first_dict_entry = tsv_dict_1[list(tsv_dict_1.keys())[0]]
    print("first_dict_entry: " + str(first_dict_entry))

    return tsv_dict_1

#this assumes that the first column of each header is the same primary key
def combine_headers(headers_1, headers_2):
    print("Combining :" + str(headers_1))
    print("With :" + str(headers_2))
    if headers_1[0] == headers_2[0]:
        exit("First column of the headers do not match. Aborting")

    combined_headers = []
    headers_2_without_primary_key = headers_2[1:]

    combined_headers = headers_1 + headers_2_without_primary_key
    print("combined_headers list: " + str(combined_headers))
    return combined_headers

def dict_to_tsv(combined_dict, combined_headers, output_tsv_path):
    with open(output_tsv_path, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(combined_headers)
        for key in combined_dict:
            row_data = [key] + combined_dict[key]
            tsv_writer.writerow(row_data)


# Main Function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('citation_tsv_path', metavar='citation_tsv_path', type=str, help='citation tsv file path')
    parser.add_argument('wiki_tsv_path', metavar='wiki_tsv_path', type=str, help='wiki tsv file path')
    parser.add_argument('output_tsv_path', metavar='output_tsv_path', type=str, help='tsv path for output file')
    args = parser.parse_args()

    citation_data, citation_headers = tsv_to_dict(args.citation_tsv_path)
    wiki_data, wiki_headers = tsv_to_dict(args.wiki_tsv_path)

    combined_dict = combine_dicts(citation_data, wiki_data) # citation data has more rows than the wiki data so use that as the "base" tsv
    combined_headers = combine_headers(citation_headers, wiki_headers)

    dict_to_tsv(combined_dict, combined_headers, args.output_tsv_path)
