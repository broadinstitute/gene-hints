'''Queries Entrez gene for citation lists for each day. Adapted from : https://github.com/pkerpedjiev/gene-citation-counts.'''

import os
import re
import sys
import urllib.request as ur
import time
from tqdm import tqdm
from optparse import OptionParser
from datetime import date, timedelta, datetime

def main():
    """
    Get the gene for citation lists for each day
    """
    usage = """
    python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --start-date 2021/07/07 --end-date 2021/07/15 --output-dir creating_citation_counts_tsv/data/tmp
    """
    num_required_args=0 # TODO: make options into args
    parser = OptionParser(usage=usage)

    # Get today's date
    today = date.today().strftime("%Y/%m/%d")

    # Set the time frame
    parser.add_option('-s', '--start-date', dest='start_date', default=today, help='The lower end of the date range to search (YYYY/MM/DD)', type='str')
    parser.add_option('-e', '--end-date', dest='end_date', default=today, help='The upper range of the date range to search (YYYY/MM/DD)', type='str')

    # Set File location for the file
    parser.add_option('-o', '--output-dir', dest='output_dir', default='creating_citation_counts_tsv/data/tmp', help='Directory to output files')

    # Get all the arguments (Start time, End time, and Output file)
    (options, args) = parser.parse_args()

    # Format start and end date
    start_date = datetime.strptime(options.start_date, '%Y/%m/%d')
    end_date = datetime.strptime(options.end_date, '%Y/%m/%d')

    days = abs((start_date - end_date).days)
    print("Get data from the last " + str(abs(days)) + " days")
    print("Output will go into: " + str(options.output_dir))
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)


    # Check for errors with setting the date and output file
    if len(args) < num_required_args:
        parser.print_help()
        sys.exit(1)

    # Get start date
    date_itr = start_date

    # Go through each day
    for i in tqdm(range(days)):

        date_itr_str = datetime.strftime(date_itr, '%Y/%m/%d')
        date_itr += timedelta(days=1)

        # Create the citation lists per each day
        # Call the entrez API
        output_file = "{}/{}.tsv".format(options.output_dir, datetime.strftime(date_itr, '%Y_%m_%d'))
        print("Output to: " + output_file)
        link = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&mindate={0}&maxdate={0}&retmax=100000".format(date_itr_str)
        time.sleep(1)
        fin = ur.urlopen(link)

        # Parse the result
        text = fin.read().decode('utf-8')
        all_pmids = re.finditer(r"<Id>(?P<pmid>[0-9]+)</Id>", text)
        out_str = ""
        for pmid_match in all_pmids:
            out_str += date_itr_str.replace('/','-') + "\t" + pmid_match.group('pmid') + "\n"

        # Save the file if there were gene publications
        if len(out_str) > 0:
            with open(output_file, 'w') as fout:
                fout.write(out_str)

# Main function
if __name__ == '__main__':
    main()

