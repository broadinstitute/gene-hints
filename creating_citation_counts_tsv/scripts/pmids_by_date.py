'''Queries Entrez gene for citation lists for each day. Adapted from : https://github.com/pkerpedjiev/gene-citation-counts.'''

import os.path as op
import re
import sys
import urllib.request as ur
import time
from tqdm import tqdm
from optparse import OptionParser
from datetime import date, timedelta, datetime

def main():
    """
    Getting the gene for citation lists for each day
    """
    usage = """
    python pmids_by_date.py 
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    # Getting today date
    today = date.today().strftime("%Y/%m/%d")

    # Setting the 5 months time frame
    days = 155

    # Getting the date 5 months ago
    month_ago = (datetime.today() - timedelta(days=155)).strftime("%Y/%m/%d")

    # Setting the time frame from today to 5 months ago
    parser.add_option('-s', '--startdate', dest='start_date', default=month_ago, help='The lower end of the date range to search (YYYY/MM/DD)', type='str')
    parser.add_option('-e', '--enddate', dest='end_date', default=today, help='The upper range of the date range to search (YYYY/MM/DD)', type='str')
    
    # Setting File location for the file
    parser.add_option('-o', '--output-dir', dest='output_dir', default='./data/pmid_by_date', help='The directory to dump all the files')

    # Getting all the arguments (Start time, End time, and Output file)
    (options, args) = parser.parse_args()

    # Formatting start and end date
    start_date = datetime.strptime(options.start_date, '%Y/%m/%d')
    end_date = datetime.strptime(options.end_date, '%Y/%m/%d')

    # Checking for errors with setting the date and output file 
    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    # Getting the start date    
    curr_date = start_date
    
    # Going through each day
    for i in tqdm(range(days)):
        # Going through each date for 5 months
        curr_date_str = datetime.strftime(curr_date, '%Y/%m/%d')
        curr_date += timedelta(days=1)

        # Creating the citation lists per each day
        # Calling the entrez API
        output_file = "creating_citation_counts_tsv/data/tmp-ssv/{}.ssv".format(datetime.strftime(curr_date, '%Y_%m_%d'))
        link = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&mindate={0}&maxdate={0}&retmax=100000".format(curr_date_str)
        time.sleep(1)
        fin = ur.urlopen(link)

        # Parsing the result
        text = fin.read().decode('utf-8')
        all_pmids = re.finditer(r"<Id>(?P<pmid>[0-9]+)</Id>", text)
        out_str = ""
        for pmid_match in all_pmids:
            out_str += curr_date_str.replace('/','-') + " " + pmid_match.group('pmid') + "\n"

        # Saving the file if there was no gene publications on a day
        if len(out_str) > 0:
            with open(output_file, 'w') as fout:
                fout.write(out_str)

# Main function
if __name__ == '__main__':
    main()

