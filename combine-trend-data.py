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
    with open(tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        isHeaderRow = True
        for row in rd:
            if isHeaderRow: # skip header
                isHeaderRow = False
                continue
            gene_symbol = row[0]
            remaining_rows = row[1:] # everything after the gene_symbol
            tsv_dict[gene_symbol] = remaining_rows

    print("row count: " + str(len(tsv_dict.keys())))
    first_dict_entry = tsv_dict[list(tsv_dict.keys())[0]]
    print("first_dict_entry: " + str(first_dict_entry))

    return tsv_dict

def combine_dicts(tsv_dict_1, tsv_dict_2):
    print("starting combine_dicts")

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

def dict_to_tsv(combined_dict, output_tsv_path):
    with open(output_tsv_path, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        headers = ["gene_symbol", "chromosome", "start", "length", "color", "full_name", "days_in_timeframe", "recent_timeframe_citation_count", "past_timeframe_citation_count", "citation_count_delta", "significance", "wikipedia_daily_page_views", "wikipedia_daily_page_views_change_from_previous_day"]
        print("combined_tsv headers list: " + str(headers))
        tsv_writer.writerow(headers)
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

    citation_data = tsv_to_dict(args.citation_tsv_path)
    wiki_data = tsv_to_dict(args.wiki_tsv_path)

    combined_dict = combine_dicts(citation_data, wiki_data) # citation data has more rows than the wiki data so use that as the "base" tsv

    dict_to_tsv(combined_dict, args.output_tsv_path)
