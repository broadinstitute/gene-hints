""" 
Creating TSVs containing the most-cited genes per species and their gene information.
"""
usage = """
    python3 creating_citation_counts_tsv/scripts/summarize_gene_citations_all_species.py ${recent_timeframe__year_pmid__ssv_path} ${past_timeframe__year_pmid__ssv_path} $days_in_timeframe
    """
    
import os
import csv
import argparse

def get_recent_gene_citation_count(species_gene2pubmed_tsv, recent_citations_ssv):
    """ 
    Input:
        species_gene2pubmed_tsv: tsv of genes and associated publication citation ID for a species
        recent_citations_ssv: ssv of all recent citations
    Output:
        Count of recent publication citations for each gene ID.
    """

    # Create a dict mapping pubmed_id to gene_id
    pubmed_gene_dict = {}
    # File header and example data row:
        # #tax_id GeneID  PubMed_ID
        # 9606    9       9173883
    with open(species_gene2pubmed_tsv) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            if row[0][0] == '#':
                continue
            gene_id = row[1]
            pubmed_id = row[2]
            if pubmed_id in pubmed_gene_dict.keys():
                gene_set = pubmed_gene_dict[pubmed_id]
                gene_set.add(gene_id)
            else:
                pubmed_gene_dict[pubmed_id] = set([gene_id])

    # Outputing the first element in pubmed_gene_dict
    print( next(iter( pubmed_gene_dict.items() )) )
    
    # Figure out which pubmed_ids in pubmed_gene_dict were published recently
    recent_pubmed_array = []
    # This file has no headers so here's an example instead:
    # 2021 34185482
    with open(recent_citations_ssv) as fd:
        rd = csv.reader(fd, delimiter=' ')
        for row in rd:
            if row[0][0] == '#':
                continue
            pubmed_id = str(row[1])
            if pubmed_id in pubmed_gene_dict.keys():
                recent_pubmed_array.append((pubmed_id))
    
    # Outputing the first element in the recent_pubmed_array variable
    if len(recent_pubmed_array) > 0:
        print(recent_pubmed_array[0])
    else:
        print("No pubmed matches -- try increasing the date range.") # todo: refactor - handle empty array gracefully

    # reverse pubmed_gene_dict so that we have a mapping of gene_id to pubmed_id
    recent_gene_pubmed_dict = {}
    for pubmed_id in recent_pubmed_array:
        gene_set = pubmed_gene_dict[pubmed_id]
        for gene_id in gene_set:
            if gene_id in recent_gene_pubmed_dict.keys():
                pubmed_set = recent_gene_pubmed_dict[gene_id]
                pubmed_set.add(pubmed_id)
            else:
                recent_gene_pubmed_dict[gene_id] = set(pubmed_id)
    
    # Dict mapping gene_id to pubmed_id citation count
    gene_citation_counts = {}
    for gene_id in recent_gene_pubmed_dict:
        gene_citation_counts[gene_id] = len(recent_gene_pubmed_dict[gene_id])

    # Outputing the first element in the gene_citation_counts variable
    print( next(iter( gene_citation_counts.items() )) )
    
    # Returning the gene_citation_counts variable
    return gene_citation_counts

# ranks a particular key in the dict (in the dict)
# ties are accounted for so that if there are 7 ties for the highest count, then all 7 genes are rank #1. The next highest gene gets rank #8.
def rank_citation_counts(gene_citation_counts):
    #first get `count`s and number of times that `count` appears
    count__num_occurances__dict = dict([])
    for key_itr in gene_citation_counts:
        count = gene_citation_counts[key_itr]
        count__num_occurances__dict[count] = count__num_occurances__dict.get(count, 0) + 1
        
    print(count__num_occurances__dict)
    sorted_count_list = sorted(count__num_occurances__dict, reverse=True)
    print("sorted_count_list: " + str(sorted_count_list))

    #calculate ranks, accounting for ties
    count__rank__dict = dict([])
    rank_counter = 1
    for count in sorted_count_list:
        count__rank__dict[count] = rank_counter
        rank_counter = rank_counter + count__num_occurances__dict[count]
    print("citation_count mapping to rank: " + str(count__rank__dict))

    #put it into a new dict for simplicity
    gene_rank_dict = {}
    for key_itr in gene_citation_counts:
        count = gene_citation_counts[key_itr]
        gene_rank_dict[key_itr] = count__rank__dict[count]
    
    return gene_rank_dict

def get_gene_info(species_gene_info_tsv, recent_gene_citation_counts, past_gene_citation_counts, recent_gene_ranks, past_gene_ranks, taxonomy_name_tsv):
    """ 
    Getting a list of the gene information per gene and species' taxonomy name.
    """

    # necessary for calculating rank delta correctly for unranked items
    recent_last_rank = int(max(recent_gene_ranks.values())) + 1
    past_last_rank = int(max(past_gene_ranks.values())) + 1
    
    # Adding gene's information to the list created from get_recent_gene_citation_count per gene and saving it to the gene_info variable
    gene_symbol__gene_info__dict = {}
    # File header and example data row:
        # #tax_id GeneID  Symbol  LocusTag        Synonyms        dbXrefs chromosome      map_location    description     type_of_gene    Symbol_from_nomenclature_authority      Full_name_from_nomenclature_authority   Nomenclature_status     Other_designations      Modification_date       Feature_type
        # 9606    59272   ACE2    -       ACEH    MIM:300335|HGNC:HGNC:13557|Ensembl:ENSG00000130234      X       Xp22.2  angiotensin converting enzyme 2 protein-coding  ACE2    angiotensin converting enzyme 2 O       angiotensin-converting enzyme 2|ACE-related carboxypeptidase|angiotensin I converting enzyme (peptidyl-dipeptidase A) 2|angiotensin I converting enzyme 2|angiotensin-converting enzyme homolog|angiotensin-converting enzyme-related carboxypeptidase|metalloprotease MPROT15|peptidyl-dipeptidase A|truncated angiotensin converting enzyme 2 20210711        -
    with open(species_gene_info_tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            if row[0][0] == '#':
                continue
            gene_id = str(row[1])
            tax_id = str(row[0])
            symbol = str(row[2])
            chromosome = str(row[6])
            full_name = str(row[8])
            recent_citation_count = recent_gene_citation_counts.get(gene_id, 0)
            past_citation_count = past_gene_citation_counts.get(gene_id, 0)
            citation_count_delta = recent_citation_count - past_citation_count
            recent_gene_rank = recent_gene_ranks.get(gene_id, recent_last_rank)
            past_gene_rank = past_gene_ranks.get(gene_id, past_last_rank)
            rank_delta = recent_gene_rank - past_gene_rank
            gene_symbol__gene_info__dict[symbol] = {
                "symbol": symbol,
                "gene_id": gene_id,
                "chromosome": chromosome,
                "full_name": full_name,
                "recent_citation_count": recent_citation_count,
                "past_citation_count": past_citation_count,
                "citation_count_delta": citation_count_delta,
                "tax_id": tax_id,
                "recent_gene_rank": recent_gene_rank,
                "past_gene_rank": past_gene_rank,
                "rank_delta": rank_delta
            }
    
    # Outputing the first element in the gene_info variable
    first_dict_entry = gene_symbol__gene_info__dict[list(gene_symbol__gene_info__dict.keys())[0]]
    print(first_dict_entry)
    
    # Getting the taxonomy ID for the species and saving it to the species_tax_id variable
    species_tax_id = first_dict_entry["tax_id"]

    # Outputing the species_tax_id variable
    print("species_tax_id: " + species_tax_id)

    # Getting the taxonomy scientific name for the taxonomy ID given and saving it to the taxonomy_scientific_name variable
    taxonomy_scientific_name = ""
    # This file has no headers so here's an example instead:
    #9606    |       Homo sapiens    |               |       scientific name |
    with open(taxonomy_name_tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            row_tax_id = str(row[0])
            name = str(row[2]) # the file does not include column headers so these are guesses based on the data
            name_type = str(row[6])
            if row_tax_id==species_tax_id and name_type=="scientific name":
                taxonomy_scientific_name = name
                break

    # Outputing the taxonomy_scientific_name variable
    print("taxonomy_scientific_name: " + taxonomy_scientific_name)

    # Returning the gene_info and taxonomy_scientific_name variable
    return gene_symbol__gene_info__dict, taxonomy_scientific_name

def get_ref_gene(species, gene_symbol__gene_info__dict): 
    """ 
    Getting a list of the gene reference information per gene.
    """

    def extract_value_from_unparsed_string(unparsed_string, key):
        """
        Extracts the value from an unparsed string
        
        key example: "gene_id"
        """
        start_index = unparsed_string.find(key) + len(key + ' "')
        end_index = unparsed_string[start_index:].find('"') + start_index
        return unparsed_string[start_index:end_index]

    
    # Going through all the files in the species directory 
    ref_gene = {}
    for file in os.listdir(f"creating_citation_counts_tsv/data/{species}"):
        #Since the reference file is uniquely named and could not be hard-coded, this will filter out files that are not .gtf files to read the reference file
        if file.endswith(".gtf"):
            # Adding gene's reference information to the list created from get_gene_info per gene and saving it to the refGene variable
            # This file has no headers so here's an example instead:
            #chr6    refGene transcript      26086290        26091034        .       -       .       gene_id "LOC108783645"; transcript_id "NR_144383";  gene_name "LOC108783645";
            # which parses to:
            #['chr6', 'refGene', 'transcript', '26086290', '26091034', '.', '-', '.', 'gene_id "LOC108783645"; transcript_id "NR_144383";  gene_name "LOC108783645";']
            with open(f"creating_citation_counts_tsv/data/{species}/{file}") as fd:
                rd = csv.reader(fd, delimiter='\t')
                for row in rd:
                    type = str(row[2]) # the file does not include column headers so these are guesses based on the data (for example: transcript, exon, 3utr, cds ...)
                    details = str(row[8])
                    if type=="transcript":
                        gene_symbol = extract_value_from_unparsed_string(details, "gene_id")
                        start_coordinate = int(row[3])
                        end_coordinate = int(row[4])
                        coordinate_length = abs(end_coordinate - start_coordinate)
                        color = "#73af42"
                        if gene_symbol in gene_symbol__gene_info__dict.keys():
                            chromosome = gene_symbol__gene_info__dict[gene_symbol]["chromosome"]
                            full_name = gene_symbol__gene_info__dict[gene_symbol]["full_name"]
                            recent_citation_count = gene_symbol__gene_info__dict[gene_symbol]["recent_citation_count"]
                            past_citation_count = gene_symbol__gene_info__dict[gene_symbol]["past_citation_count"]
                            citation_count_delta = gene_symbol__gene_info__dict[gene_symbol]["citation_count_delta"]
                            recent_gene_rank = gene_symbol__gene_info__dict[gene_symbol]["recent_gene_rank"]
                            past_gene_rank = gene_symbol__gene_info__dict[gene_symbol]["past_gene_rank"]
                            rank_delta = gene_symbol__gene_info__dict[gene_symbol]["rank_delta"]
                        else:
                            chromosome = ""
                            full_name = ""
                            recent_citation_count = 0
                            past_citation_count = 0
                            citation_count_delta = 0
                            recent_gene_rank = -1
                            past_gene_rank = -1
                            rank_delta = -1
                        ref_gene[gene_symbol] = {
                            "chromosome": chromosome,
                            "start_coordinate": start_coordinate,
                            "coordinate_length": coordinate_length,
                            "color": color,
                            "full_name": full_name,
                            "recent_citation_count": recent_citation_count,
                            "past_citation_count": past_citation_count,
                            "citation_count_delta": citation_count_delta,
                            "recent_gene_rank": recent_gene_rank,
                            "past_gene_rank": past_gene_rank,
                            "rank_delta": rank_delta
                            }

    # Outputing the ref_gene variable
    print( next(iter( ref_gene.items() )) )
    
    # Returning the ref_gene variable
    return ref_gene

def sort_and_list_genes(ref_gene):
    """
    Returns a list of the genes in order of citation count, starting with the largest first.
    """

    # Getting the list of the top ten most cited genes from the list created from get_ref_gene
    sorted_genes_list = sorted(ref_gene.items(), key=lambda x: x[1]['citation_count_delta'], reverse=True)
    
    print("Gene with the greatest positive citation delta: " + str(sorted_genes_list[0]))

    return sorted_genes_list

def create_tsv_for_genes(sorted_genes_list, tax_name, timeframe_days):
    """
    Creating TSVs containing gene information.
    """
    
    # Getting the species' scientific taxonomy name and saving it to the gene_species_name variable
    gene_species_name = tax_name.replace(" ", "-").replace("(", "-").replace(")", "-").replace("/", "-").replace("=", "-").lower()
    
    # Outputting the gene_species_name variable
    print(gene_species_name)
    
    # Creating the TSV name that contains the gene_species_name variable
    tsv_name= f'data/{gene_species_name}-citation-information.tsv'
    
    # Creating the TSV
    with open(tsv_name, 'wt') as out_file:
        # Setting up the TSV
        tsv_writer = csv.writer(out_file, delimiter='\t')
        
        # Add the header row to the TSV
        tsv_writer.writerow(["gene_symbol", "chromosome", "start", "length", "color", "full_name", "days_in_timeframe", "recent_timeframe_citation_count", "past_timeframe_citation_count", "citation_count_delta", "significance", "recent_gene_rank", "past_gene_rank", "rank_delta",])
        
        # Going through each gene in the sorted_genes_list 
        for gene in sorted_genes_list:
            # Adding the gene information to the TSV
            symbol = gene[0]
            chromosome = gene[1]["chromosome"]
            start = gene[1]["start_coordinate"]
            length = gene[1]["coordinate_length"]
            color = gene[1]["color"]
            full_name = gene[1]["full_name"]
            recent_citation_count = gene[1]["recent_citation_count"]
            past_citation_count = gene[1]["past_citation_count"]
            citation_count_delta = gene[1]["citation_count_delta"]
            significance = "" # `significance` is now being handled from the front-end. leaving this here for backwards compatibility. Remove when no longer needed.
            recent_gene_rank = gene[1]["recent_gene_rank"]
            past_gene_rank = gene[1]["past_gene_rank"]
            rank_delta = gene[1]["rank_delta"]
            tsv_row_values = [symbol, chromosome, start, length, color, full_name, timeframe_days, recent_citation_count, past_citation_count, citation_count_delta, significance, recent_gene_rank, past_gene_rank, rank_delta]

            tsv_writer.writerow(tsv_row_values)
    
# Main Function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('recent_citation_count_ssv', metavar='recent_citation_count_ssv', type=str, help='recent citation count ssv file path')
    parser.add_argument('past_citation_count_ssv', metavar='past_citation_count_ssv', type=str, help='past citation count ssv file path')
    parser.add_argument('timeframe_days', metavar='timeframe_days', type=int, help='days in the timeframe')
    args = parser.parse_args()
    
    # Setting list of species
    list_of_species = ["human", "mouse", "rat", "dog", "cat"]

    # Going through each species
    for species in list_of_species:
        # Outputing species name
        print(species)
        
        recent_gene_citation_counts = get_recent_gene_citation_count(f"creating_citation_counts_tsv/data/{species}/gene2pubmed", args.recent_citation_count_ssv)
        past_gene_citation_counts = get_recent_gene_citation_count(f"creating_citation_counts_tsv/data/{species}/gene2pubmed", args.past_citation_count_ssv)

        recent_gene_rank = rank_citation_counts(recent_gene_citation_counts)
        past_gene_rank = rank_citation_counts(past_gene_citation_counts)
        
        gene_symbol__gene_info__dict, tax_name = get_gene_info(f"creating_citation_counts_tsv/data/{species}/gene_info", recent_gene_citation_counts, past_gene_citation_counts, recent_gene_rank, past_gene_rank, "creating_citation_counts_tsv/taxonomy_name")
        
        ref_gene = get_ref_gene(species, gene_symbol__gene_info__dict)

        sorted_genes_list = sort_and_list_genes(ref_gene)

        create_tsv_for_genes(sorted_genes_list, tax_name, args.timeframe_days)
