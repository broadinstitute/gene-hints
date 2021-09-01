"""
Create TSVs containing the most-cited genes per species and their gene information.
"""
usage = """
    python3 creating_citation_counts_tsv/scripts/summarize_gene_citations_all_species.py ${pmid_times_path} ${prev_pmid_times_path} $days_in_timeframe
    """

import os
import csv
import argparse

organisms = [
    ["homo-sapiens", "9606", "hg38"], # human
    ["mus-musculus", "10090", "mm39"], # mouse
    ["rattus-norvegicus", "10116", "rn6"], # rat
    ["canis-lupus-familiaris", "9615", "canFam5"], # dog
    ["felis-catus", "9685", "felCat9"], # cat
]

cites_dir = './pubmed_citations/'
data_dir = cites_dir + 'data/'

def get_cites_by_gene(gene2pubmed_tsv, cites_tsv):
    """
    Input:
        species_gene2pubmed_tsv: tsv of genes and associated publication cite ID for a species
        cites_tsv: tsv of all recent citations
    Output:
        Count of recent publication citations for each gene ID.
    """

    # Map PubMed IDs (PMIDs) to gene IDs
    genes_by_pmid  = {}
    # File header and example data row:
        # #tax_id GeneID  pmid
        # 9606    9       9173883
    print(gene2pubmed_tsv)
    with open(gene2pubmed_tsv) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            if row[0][0] == '#':
                continue
            gene_id = row[1]
            pmid = row[2]
            if pmid in genes_by_pmid.keys():
                gene_set = genes_by_pmid[pmid]
                gene_set.add(gene_id)
            else:
                genes_by_pmid[pmid] = set([gene_id])

    # # Output the first element in genes_by_pmid
    print('genes_by_pmid')
    print(next(iter( genes_by_pmid.items() )) )

    # Figure out which PMIDs in genes_by_pmid were published recently
    pmids = []
    # This file has no headers so here's an example instead:
    # 2021 34185482
    with open(cites_tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            if row[0][0] == '#' or len(row) < 2:
                print('Short row' + str(row))
                continue
            pmid = str(row[1])
            if pmid in genes_by_pmid.keys():
                pmids.append((pmid))

    # Output the first element in the pmids variable
    if len(pmids) > 0:
        print('pmids[0]')
        print(pmids[0])
    else:
        print("No PubMed matches -- try increasing the date range.") # todo: refactor - handle empty array gracefully

    # Reverse genes_by_pmid so that we have a mapping of gene_id to pmid
    pmids_by_gene = {}
    for pmid in pmids:
        gene_set = genes_by_pmid[pmid]
        for gene_id in gene_set:
            if gene_id in pmids_by_gene.keys():
                pubmed_set = pmids_by_gene[gene_id]
                pubmed_set.add(pmid)
            else:
                pmids_by_gene[gene_id] = set(pmid)

    # Dict mapping gene_id to PMID citation count
    cites_by_gene = {}
    for gene_id in pmids_by_gene:
        cites_by_gene[gene_id] = len(pmids_by_gene[gene_id])

    # Output the first element in the cites variable
    print('len(cites_by_gene.items())')
    print(len(cites_by_gene.items()))
    print('cites[0]')
    print( next(iter( cites_by_gene.items() )) )

    return cites_by_gene

def rank_counts(counts_by_key):
    """
    Ranks a particular key in the dict

    For ties, if there are 7 ties for the highest count, then all 7 keys
    (e.g. genes) are rank #1. The next highest key gets rank #8.
    """
    # First get `count`s and number of times that `count` appears
    occurences_by_count = dict([])
    for key in counts_by_key:
        count = counts_by_key[key]
        occurences_by_count[count] = occurences_by_count.get(count, 0) + 1

    sorted_occurences = sorted(occurences_by_count, reverse=True)
    print("sorted occurences_by_count: " + str(sorted_occurences))

    # Calculate ranks, accounting for ties
    ranks_by_count = dict([])
    rank_counter = 1
    for count in sorted_occurences:
        ranks_by_count[count] = rank_counter
        rank_counter = rank_counter + occurences_by_count[count]
    print("cite mapping to rank: " + str(ranks_by_count))

    # Put it into a new dict for simplicity
    ranks_by_key = {}
    for key in counts_by_key:
        count = counts_by_key[key]
        ranks_by_key[key] = ranks_by_count[count]

    return ranks_by_key

def enrich_gene_info(gene_info_file, cites_by_gene, prev_cites_by_gene, cite_ranks, prev_cite_ranks):
    """
    Gets a list of the gene information per gene
    """

    # Needed to calculate rank delta for unranked items
    last_rank = int(max(cite_ranks.values())) + 1
    prev_last_rank = int(max(prev_cite_ranks.values())) + 1

    # Add gene's information to the list created from get_cite
    # per gene and save it to the gene_info variable
    gene_info_by_symbol = {}
    # File header and example data row:
        # #tax_id GeneID  Symbol  LocusTag        Synonyms        dbXrefs chromosome      map_location    description     type_of_gene    Symbol_from_nomenclature_authority      Full_name_from_nomenclature_authority   Nomenclature_status     Other_designations      Modification_date       Feature_type
        # 9606    59272   ACE2    -       ACEH    MIM:300335|HGNC:HGNC:13557|Ensembl:ENSG00000130234      X       Xp22.2  angiotensin converting enzyme 2 protein-coding  ACE2    angiotensin converting enzyme 2 O       angiotensin-converting enzyme 2|ACE-related carboxypeptidase|angiotensin I converting enzyme (peptidyl-dipeptidase A) 2|angiotensin I converting enzyme 2|angiotensin-converting enzyme homolog|angiotensin-converting enzyme-related carboxypeptidase|metalloprotease MPROT15|peptidyl-dipeptidase A|truncated angiotensin converting enzyme 2 20210711        -
    with open(gene_info_file) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            if row[0][0] == '#':
                continue

            gene_id = str(row[1])
            taxid = str(row[0])
            symbol = str(row[2])
            chromosome = str(row[6])
            full_name = str(row[8])
            cites = cites_by_gene.get(gene_id, 0)
            prev_cites = prev_cites_by_gene.get(gene_id, 0)
            cite_delta = cites - prev_cites
            cite_rank = cite_ranks.get(gene_id, last_rank)
            prev_cite_rank = prev_cite_ranks.get(gene_id, prev_last_rank)
            rank_delta = cite_rank - prev_cite_rank

            gene_info_by_symbol[symbol] = {
                "symbol": symbol,
                "gene_id": gene_id,
                "chromosome": chromosome,
                "full_name": full_name,
                "cites": cites,
                "prev_cites": prev_cites,
                "cite_delta": cite_delta,
                "taxid": taxid,
                "cite_rank": cite_rank,
                "prev_cite_rank": prev_cite_rank,
                "rank_delta": rank_delta
            }

    return gene_info_by_symbol

def extract_value_from_unparsed_string(unparsed_string, key):
    """
    Extracts the value from an unparsed string

    key example: "gene_id"
    """
    start_index = unparsed_string.find(key) + len(key + ' "')
    end_index = unparsed_string[start_index:].find('"') + start_index
    return unparsed_string[start_index:end_index]

def get_ref_gene(organism, gene_info_by_symbol):
    """
    Gets a list of the gene reference information per gene.
    """

    # Iterate files in the species directory
    ref_gene = {}
    print('data_dir + organism', data_dir + organism)
    for file in os.listdir(data_dir + organism):
        # Since the reference file is uniquely named and could not be
        # hard-coded, this will filter out files that are not .gtf files to
        # read the reference file
        if not file.endswith(".gtf"):
            continue

        # Add gene's reference information to the list created from
        # enrich_gene_info per gene and saving it to the refGene variable
        # This file has no headers so here's an example instead:
        #chr6    refGene transcript      26086290        26091034        .       -       .       gene_id "LOC108783645"; transcript_id "NR_144383";  gene_name "LOC108783645";
        # which parses to:
        #['chr6', 'refGene', 'transcript', '26086290', '26091034', '.', '-', '.', 'gene_id "LOC108783645"; transcript_id "NR_144383";  gene_name "LOC108783645";']
        org_geneinfo = f"{data_dir}{organism}/{file}"
        with open(org_geneinfo) as fd:
            rd = csv.reader(fd, delimiter='\t')
            for row in rd:
                # File lacks column headers so these are guesses based on the data
                # E.g.: transcript, exon, 3utr, cds ...
                type = str(row[2])
                details = str(row[8])
                if type != "transcript":
                    continue

                gene_symbol = extract_value_from_unparsed_string(details, "gene_id")
                start_coordinate = int(row[3])
                end_coordinate = int(row[4])
                coordinate_length = abs(end_coordinate - start_coordinate)
                color = "#73af42"

                if gene_symbol in gene_info_by_symbol.keys():
                    gene_info = gene_info_by_symbol[gene_symbol]
                    chromosome = gene_info["chromosome"]
                    full_name = gene_info["full_name"]
                    cites = gene_info["cites"]
                    prev_cites = gene_info["prev_cites"]
                    cite_delta = gene_info["cite_delta"]
                    cite_rank = gene_info["cite_rank"]
                    prev_cite_rank = gene_info["prev_cite_rank"]
                    rank_delta = gene_info["rank_delta"]
                else:
                    chromosome = ""
                    full_name = ""
                    cites = 0
                    prev_cites = 0
                    cite_delta = 0
                    cite_rank = -1
                    prev_cite_rank = -1
                    rank_delta = -1

                ref_gene[gene_symbol] = {
                    "chromosome": chromosome,
                    "start_coordinate": start_coordinate,
                    "coordinate_length": coordinate_length,
                    "color": color,
                    "full_name": full_name,
                    "cites": cites,
                    "prev_cites": prev_cites,
                    "cite_delta": cite_delta,
                    "cite_rank": cite_rank,
                    "prev_cite_rank": prev_cite_rank,
                    "rank_delta": rank_delta
                }

    # Print ref_gene variable
    print('ref_gene')
    print( next(iter( ref_gene.items() )) )

    return ref_gene

def sort_and_list_genes(ref_gene):
    """
    Returns a list of the genes in order of cite count, starting with the largest first
    """

    # Get the list of the top ten most cited genes from the list created from get_ref_gene
    sorted_genes_list = sorted(ref_gene.items(), key=lambda x: x[1]['cite_delta'], reverse=True)

    print("Gene with the greatest positive cite delta: " + str(sorted_genes_list[0]))

    return sorted_genes_list

def create_tsv_for_genes(sorted_genes_list, org_name, timeframe_days):
    """
    Creates TSVs containing gene information
    """

    # Simplify organism name to ease use as file name
    safe_org_name = org_name.replace(" ", "-").replace("(", "-").replace(")", "-").replace("/", "-").replace("=", "-").lower()

    print(safe_org_name)

    # Create the TSV name that contains the gene_species_name variable
    tsv_name= f'data/{safe_org_name}-citation-information.tsv'

    # Create the TSV
    with open(tsv_name, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')

        # Add the header row to the TSV
        tsv_writer.writerow([
            "gene_symbol", "chromosome", "start", "length", "color",
            "full_name", "days_in_timeframe",
            "cites",
            "prev_cites",
            "citation_delta",
            "cite_rank",
            "prev_cite_rank",
            "rank_delta"])

        # Go through each gene in the sorted_genes_list
        for gene in sorted_genes_list:
            # Add the gene information to the TSV
            symbol = gene[0]
            chromosome = gene[1]["chromosome"]
            start = gene[1]["start_coordinate"]
            length = gene[1]["coordinate_length"]
            color = gene[1]["color"]
            full_name = gene[1]["full_name"]
            cites = gene[1]["cites"]
            prev_cites = gene[1]["prev_cites"]
            cite_delta = gene[1]["cite_delta"]
            cite_rank = gene[1]["cite_rank"]
            prev_cite_rank = gene[1]["prev_cite_rank"]
            rank_delta = gene[1]["rank_delta"]
            tsv_row_values = [
                symbol, chromosome, start, length, color, full_name,
                timeframe_days, cites, prev_cites,
                cite_delta, cite_rank, prev_cite_rank,
                rank_delta
            ]

            tsv_writer.writerow(tsv_row_values)

# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument(
        'citation_tsv',
        metavar='citation_tsv',
        type=str,
        help='Path to file containing citation counts over timeframe'
    )
    parser.add_argument(
        'prev_citation_tsv',
        metavar='prev_citation_tsv',
        type=str,
        help='Path to file containing citation counts over previous timeframe'
    )
    parser.add_argument(
        'timeframe_days',
        metavar='timeframe_days',
        type=int,
        help='Days in the timeframe'
    )
    args = parser.parse_args()

    for org_array in organisms:
        organism = org_array[0]
        print(organism)

        org_dir = data_dir + organism

        cites_by_gene = get_cites_by_gene(f"{org_dir}/gene2pubmed", args.citation_tsv)
        prev_cites_by_gene = get_cites_by_gene(
            f"{org_dir}/gene2pubmed",
            args.prev_citation_tsv
        )

        cite_rank = rank_counts(cites_by_gene)
        prev_cite_rank = rank_counts(cites_by_gene)

        gene_info = enrich_gene_info(
            f"{org_dir}/gene_info",
            cites_by_gene,
            prev_cites_by_gene,
            cite_rank,
            prev_cite_rank
        )

        scientific_name = organism

        ref_gene = get_ref_gene(organism, gene_info)

        sorted_genes_list = sort_and_list_genes(ref_gene)

        create_tsv_for_genes(sorted_genes_list, scientific_name, args.timeframe_days)
