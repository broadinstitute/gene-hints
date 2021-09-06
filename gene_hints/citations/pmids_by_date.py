"""Query PubMed for citation lists for each day

Adapted from https://github.com/pkerpedjiev/gene-citation-counts
"""

import os
import re
import sys
import urllib.request as ur
import time
from tqdm import tqdm
from optparse import OptionParser
from datetime import date as dt_date, timedelta, datetime

def pmids_by_date(
    start_date=None, end_date=None,
    output_dir="creating_citation_counts_tsv/data/tmp"
):
    """Query PubMed for citation lists for each day
    """

    # Get today's date
    today = dt_date.today().strftime("%Y/%m/%d")

    if start_date == None:
        start_date = today
    if end_date == None:
        end_date = today

    # Format start and end date
    start_date = datetime.strptime(start_date, "%Y/%m/%d")
    end_date = datetime.strptime(end_date, "%Y/%m/%d")

    days = abs((start_date - end_date).days)
    print("Get data from the last " + str(abs(days)) + " days")
    print("Output will go into: " + output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get start date
    date_itr = start_date

    eutils_base =\
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"

    # Go through each day
    for i in tqdm(range(days)):

        date = datetime.strftime(date_itr, "%Y/%m/%d")
        date_itr += timedelta(days=1)

        # Create the citation lists per each day
        date_underscore = date.replace("/", "_")

        # Call NCBI E-utils / Entrez API
        # https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.mindate_maxdate
        link = f"{eutils_base}&mindate={date}&maxdate={date}&retmax=100000"
        time.sleep(1)
        fin = ur.urlopen(link)

        # Parse the result
        text = fin.read().decode("utf-8")
        all_pmids = re.finditer(r"<Id>(?P<pmid>[0-9]+)</Id>", text)
        date_pmid_lines = ""
        for pmid_match in all_pmids:
            date_hyphen = date.replace("/", "-")
            pmid = pmid_match.group("pmid")
            date_pmid_lines += date_hyphen + "\t" + pmid + "\n"

        # Save the file only if there were publications
        if len(date_pmid_lines) == 0:
            continue

        output_path = f"{output_dir}/{date_underscore}.tsv"
        with open(output_path, "w") as f:
            f.write(date_pmid_lines)

        print("Output to: " + output_path)

# Command-line handler
if __name__ == "__main__":
    usage = """
    python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --start-date 2021/07/07 --end-date 2021/07/15 --output-dir creating_citation_counts_tsv/data/tmp
    """
    parser = OptionParser(usage=usage)

    # Set the time frame
    parser.add_option(
        "-s", "--start-date", dest="start_date",
        help="Recent end of the date range to search (YYYY/MM/DD)"
    )
    parser.add_option(
        "-e", "--end-date", dest="end_date",
        help="Distant end the date range to search (YYYY/MM/DD)"
    )

    parser.add_option(
        "-o", "--output-dir", dest="output_dir",
        help="Directory to output files"
    )

    # Get all the arguments (Start time, End time, and Output file)
    (options, args) = parser.parse_args()

    num_required_args = 0 # TODO: make options into args

    # Check for errors with setting the date and output file
    if len(args) < num_required_args:
        parser.print_help()
        sys.exit(1)


    pmids_by_date(options.start_date, options.end_date, options.output_dir)

