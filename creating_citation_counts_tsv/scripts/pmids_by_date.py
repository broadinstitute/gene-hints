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
    Getting the gene for citation lists for each day
    """
    usage = """
    python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --startdate 2021/07/07 --enddate 2021/07/15 --output-dir creating_citation_counts_tsv/data/tmp-ssv
    """
    num_required_args=0 #todo make options into args
    parser = OptionParser(usage=usage)

    # Getting today date
    today = date.today().strftime("%Y/%m/%d")

    # Setting the time frame
    parser.add_option('-s', '--startdate', dest='start_date', default=today, help='The lower end of the date range to search (YYYY/MM/DD)', type='str')
    parser.add_option('-e', '--enddate', dest='end_date', default=today, help='The upper range of the date range to search (YYYY/MM/DD)', type='str')
    
    # Setting File location for the file
    parser.add_option('-o', '--output-dir', dest='output_dir', default='creating_citation_counts_tsv/data/tmp-ssv', help='The directory to dump all the files')

    # Getting all the arguments (Start time, End time, and Output file)
    (options, args) = parser.parse_args()

    # Formatting start and end date
    start_date = datetime.strptime(options.start_date, '%Y/%m/%d')
    end_date = datetime.strptime(options.end_date, '%Y/%m/%d')
    
    days = abs((start_date - end_date).days)
    print("Getting data from the last " + str(abs(days)) + " days")
    print("Output will go into: " + str(options.output_dir))
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)


    # Checking for errors with setting the date and output file 
    if len(args) < num_required_args:
        parser.print_help()
        sys.exit(1)

    # Getting the start date    
    date_itr = start_date

    # Going through each day
    for i in tqdm(range(days)):
        # Going through each date
        date_itr_str = datetime.strftime(date_itr, '%Y/%m/%d')
        date_itr += timedelta(days=1)

        # Creating the citation lists per each day
        # Calling the entrez API
        output_file = "{}/{}.ssv".format(options.output_dir, datetime.strftime(date_itr, '%Y_%m_%d'))
        print("Outputting to: " + output_file)
        link = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&mindate={0}&maxdate={0}&retmax=100000".format(date_itr_str)
        time.sleep(1)
        fin = ur.urlopen(link)

        # Parsing the result
        text = fin.read().decode('utf-8')
        all_pmids = re.finditer(r"<Id>(?P<pmid>[0-9]+)</Id>", text)
        out_str = ""
        for pmid_match in all_pmids:
            out_str += date_itr_str.replace('/','-') + " " + pmid_match.group('pmid') + "\n"

        # Saving the file if there were gene publications
        if len(out_str) > 0:
            with open(output_file, 'w') as fout:
                fout.write(out_str)

# Main function
if __name__ == '__main__':
    main()

