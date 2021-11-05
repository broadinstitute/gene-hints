"""Combine genomic and citation data, compute statistics, and write to TSV

Inspired by https://github.com/pkerpedjiev/gene-citation-counts
"""

import csv
import argparse

from lib import read_organisms

def get_pmids_with_genes_in_timeframe(genes_by_pmid, pmid_dates_path):
    """Return PMIDs in genes_by_pmid that were published in a timeframe

    The timeframe is implied by data in the file at `pmid_dates_path`.
    """
    pmids = []

    # This file has no headers so here's an example instead:
    # 2021 34185482
    with open(pmid_dates_path) as fd:
        reader = csv.reader(fd, delimiter="\t")
        for row in reader:
            if row[0][0] == "#":
                continue
            pmid = str(row[1])
            if pmid in genes_by_pmid:
                pmids.append((pmid))

    return pmids

def rank_counts(counts_by_key):
    """Rank keys in a dict by integer count values

    If there are 7 ties for the highest count, then all 7 keys
    (e.g. genes) are rank #1. The next highest key gets rank #8.
    """
    # First get `count`s and number of times that `count` appears
    occurences_by_count = dict([])
    for key in counts_by_key:
        count = counts_by_key[key]
        occurences_by_count[count] = occurences_by_count.get(count, 0) + 1

    sorted_occurences = sorted(occurences_by_count, reverse=True)
    # print("sorted occurences_by_count: " + str(sorted_occurences))

    # Calculate ranks, accounting for ties
    ranks_by_count = dict([])
    rank_counter = 1
    for count in sorted_occurences:
        ranks_by_count[count] = rank_counter
        rank_counter = rank_counter + occurences_by_count[count]
    # print("cite mapping to rank: " + str(ranks_by_count))

    # Put it into a new dict for simplicity
    ranks_by_key = {}
    for key in counts_by_key:
        count = counts_by_key[key]
        ranks_by_key[key] = ranks_by_count[count]

    return ranks_by_key

def sort_genes(gene_dict, key):
    """Return list of genes ordered by key, in descending order
    """
    sorted_genes = sorted(
        gene_dict.items(), key=lambda x: x[1][key],
        reverse=True
    )

    # print(f"Gene with highest {key}: {str(sorted_genes[0])}")

    return sorted_genes

def pretty_print_table(raw_rows, num_rows):
    """Print rows with left-justified columns, for easy reading
    """
    print(f"First {num_rows} rows")
    rows = raw_rows[1: num_rows] # Offset accounts for meta-information

    print(raw_rows[0][0])
    # From https://stackoverflow.com/a/12065663
    widths = [max(map(len, col)) for col in zip(*rows)]
    for row in rows:
        print("  ".join((val.ljust(width) for val, width in zip(row, widths))))

def pretty_org_name(org_name):
    """Convert e.g. "homo-sapiens" to "Homo sapiens"
    """
    first_letter = org_name[0].upper()
    return first_letter + org_name[1:].replace("-", " ")


class EnrichCitations():
    def __init__(
        self,
        data_dir="./pubmed_citations/data/",
        output_dir="./data/"
    ):
        self.data_dir = data_dir
        self.output_dir = output_dir

    def get_output_path(self, organism):
        return f"{self.output_dir}/{organism}-pubmed-citations.tsv"

    def get_genes_by_pmid(self, organism):
        """Map PubMed IDs (PMIDs) to gene IDs
        """
        genes_by_pmid  = {}

        gene2pubmed_path = self.data_dir + organism + "/gene2pubmed"

        # File header and example data row:
            # #tax_id GeneID  pmid
            # 9606    9       9173883
        # print(gene2pubmed_path)
        with open(gene2pubmed_path) as f:
            reader = csv.reader(f, delimiter="\t", quotechar='"')
            for row in reader:
                if row[0][0] == "#":
                    continue
                gene_id = row[1]
                pmid = row[2]
                if pmid in genes_by_pmid.keys():
                    gene_set = genes_by_pmid[pmid]
                    gene_set.add(gene_id)
                else:
                    genes_by_pmid[pmid] = set([gene_id])

        return genes_by_pmid

    def get_cites_by_gene(self, organism, pmid_dates_path):
        """Get citation counts for an organism's genes over time

        Input:
            organism: scientific name of species, hyphen-case (e.g. homo-sapiens)
            pmid_dates_path: location of TSV containing daily cite counts
        Output:
            Count of citations for each gene ID.
        """

        # Get a map PubMed IDs (PMIDs) to gene IDs.
        # These PMIDs represent publications that mention a gene in the organism.
        genes_by_pmid  = self.get_genes_by_pmid(organism)

        # Output the first element in genes_by_pmid
        # print('Number of genes in first PMID in genes_by_pmid')
        # print(len(list(list(genes_by_pmid.items())[0][1])))

        # Get list of above PMIDs that were published in the dates of interest
        pmids = get_pmids_with_genes_in_timeframe(genes_by_pmid, pmid_dates_path)

        # if len(pmids) > 0:
        #     # Output the first element in the pmids variable
        #     print('pmids[0]')
        #     print(pmids[0])

        # Invert genes_by_pmid, and only include PMIDs in timeframe
        pmids_by_gene = {}
        for pmid in pmids:
            gene_set = genes_by_pmid[pmid]
            for gene_id in gene_set:
                if gene_id in pmids_by_gene.keys():
                    pubmed_set = pmids_by_gene[gene_id]
                    pubmed_set.add(pmid)
                else:
                    pmids_by_gene[gene_id] = set(pmid)

        # Collapse list of PMIDs to count of PMIDs (citation counts, "cites")
        cites_by_gene = {}
        for gene_id in pmids_by_gene:
            cites_by_gene[gene_id] = len(pmids_by_gene[gene_id])

        # cites = cites_by_gene.items()
        # # Output the first element in the cites variable
        # if (len(list(cites))) > 0:
        #     print('cites[0]')
        #     print( next(iter( cites )) )

        return cites_by_gene

    def enrich_genes(
        self, organism, cites_by_gene, prev_cites_by_gene,
        cite_ranks, prev_cite_ranks
    ):
        """Combine genomic data and citation data
        """

        gene_info_path = self.data_dir + organism + "/gene_info"

        # Needed to calculate rank delta for unranked items
        last_rank = int(max(cite_ranks.values())) + 1
        prev_last_rank = int(max(prev_cite_ranks.values())) + 1

        # Add gene's information to the list created from get_cite
        # per gene and save it to the gene_info variable
        enriched_genes = {}
        # File header and example data row:
            # #tax_id GeneID  Symbol  LocusTag        Synonyms        dbXrefs chromosome      map_location    description     type_of_gene    Symbol_from_nomenclature_authority      Full_name_from_nomenclature_authority   Nomenclature_status     Other_designations      Modification_date       Feature_type
            # 9606    59272   ACE2    -       ACEH    MIM:300335|HGNC:HGNC:13557|Ensembl:ENSG00000130234      X       Xp22.2  angiotensin converting enzyme 2 protein-coding  ACE2    angiotensin converting enzyme 2 O       angiotensin-converting enzyme 2|ACE-related carboxypeptidase|angiotensin I converting enzyme (peptidyl-dipeptidase A) 2|angiotensin I converting enzyme 2|angiotensin-converting enzyme homolog|angiotensin-converting enzyme-related carboxypeptidase|metalloprotease MPROT15|peptidyl-dipeptidase A|truncated angiotensin converting enzyme 2 20210711        -
        with open(gene_info_path) as fd:
            reader = csv.reader(fd, delimiter="\t")
            for row in reader:
                if row[0][0] == "#":
                    continue

                # Genomic data, from gene_info_file
                gene_id = str(row[1])
                symbol = str(row[2])

                if (
                    gene_id not in cites_by_gene and
                    gene_id not in prev_cites_by_gene
                ):
                    # If gene is uncited, skip it
                    continue

                # Citation data
                cites = cites_by_gene.get(gene_id, 0)
                prev_cites = prev_cites_by_gene.get(gene_id, 0)
                cite_delta = cites - prev_cites
                cite_rank = cite_ranks.get(gene_id, last_rank)
                prev_cite_rank = prev_cite_ranks.get(gene_id, prev_last_rank)

                # Rank 1 is highest, so order negative rank change higher
                cite_rank_delta = prev_cite_rank - cite_rank

                enriched_genes[symbol] = {
                    "symbol": symbol,
                    "cites": cites,
                    "cite_delta": cite_delta,
                    "cite_rank": cite_rank,
                    "cite_rank_delta": cite_rank_delta
                }

        return enriched_genes

    def save_to_file(self, sorted_genes_list, organism, days):
        """Write TSV file that combines citation and genomic data
        """
        output_path = self.get_output_path(organism)
        no_coordinates = []
        rows = []

        with open(output_path, "wt") as f:
            tsv_writer = csv.writer(f, delimiter="\t")

            meta = ["## cite_days=" + str(days)]
            tsv_writer.writerow(meta)
            rows.append(meta)

            header = [
                "# gene", "cites", "cite_delta", "cite_rank",
                "cite_rank_delta"
            ]
            # Add header row to the TSV
            tsv_writer.writerow(header)
            rows.append(header)

            # Write a row in the TSV for each gene in the sorted list
            for item in sorted_genes_list:
                # Add the gene information to the TSV
                symbol = item[0]
                gene = item[1]
                cites = gene["cites"]
                cite_delta = gene["cite_delta"]
                cite_rank = gene["cite_rank"]
                cite_rank_delta = gene["cite_rank_delta"]
                row = [
                    symbol, cites,
                    cite_delta, cite_rank,
                    cite_rank_delta
                ]
                row = [str(item) for item in row]
                rows.append(row)
                tsv_writer.writerow(row)

        num_missing = len(no_coordinates)
        if num_missing > 0:
            no_coordinates = ", ".join(no_coordinates)
            print(f"No genomic coordinates for {num_missing} genes: {no_coordinates}")
        else:
            print("All cited genes in this organism had genomic coordinates!")

        print(f"Wrote { len(rows) } gene citation hints to {output_path}")
        pretty_print_table(rows, 10)

    def run(
        self, pmid_dates_path, prev_pmid_dates_path, days, sort_by="count"
    ):
        """Construct cite hints, write to TSV.  Intended for use in other modules.
        """
        organisms = read_organisms()
        sort_map = {
            "count": "cites",
            "delta": "cite_delta",
            "rank": "cite_rank",
            "rank_delta": "cite_rank_delta"
        }
        sort_key = sort_map[sort_by]

        for org in organisms:
            organism = org["scientific_name"]
            common = org["common_name"]
            pretty_org = pretty_org_name(organism)
            print("\n")
            print(f"Enriching citations for {pretty_org} ({common})")

            cites_by_gene = self.get_cites_by_gene(organism, pmid_dates_path)
            prev_cites_by_gene = self.get_cites_by_gene(organism, prev_pmid_dates_path)

            if (
                len(list(cites_by_gene)) == 0 and
                len(list(prev_cites_by_gene)) == 0
            ):
                print(
                    f"No publications in PubMed cited genes from {organism} " +
                    f"in the last {days} days.  " +
                    f"Try increasing the `days` value in `citations.py`."
                )
                continue

            cite_rank = rank_counts(cites_by_gene)
            prev_cite_rank = rank_counts(prev_cites_by_gene)

            enriched_genes = self.enrich_genes(
                organism,
                cites_by_gene,
                prev_cites_by_gene,
                cite_rank,
                prev_cite_rank
            )

            cite_hints = sort_genes(enriched_genes, sort_key)

            self.save_to_file(cite_hints, organism, days)

# Command-line handler
if __name__ == "__main__":
    usage = """
    python3 gene_hints/citations/enrich_citations.py --pmid-dates-path ./pubmed_citations/data/pmid_dates.tsv --prev-pmid-dates-path ./pubmed_citations/data/prev_pmid_dates.tsv --days 5
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=usage
    )

    parser.add_argument(
        "--pmid-dates-path",
        help="Path to file containing citation counts over time"
    )
    parser.add_argument(
        "--prev-pmid-dates-path",
        help="Path to file containing citation counts over previous timeframe"
    )
    parser.add_argument(
        "--days",
        type=int,
        help="Days in the timeframe"
    )
    parser.add_argument(
        "--sort-by",
        help="Metric by which to sort PubMed citations",
        choices=["count", "delta", "rank", "rank_delta"],
        default="count"
    )
    args = parser.parse_args()
    pmid_dates_path = args.pmid_dates_path
    prev_pmid_dates_path = args.prev_pmid_dates_path
    days = args.days
    sort_by = args.sort_by

    enrich_citations = EnrichCitations()
    enrich_citations.run(pmid_dates_path, prev_pmid_dates_path, days, sort_by)
