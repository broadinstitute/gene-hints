""" 
Creating TSVs containing the 10 most-cited genes per species and their gene information.
"""
import os
from datetime import date, timedelta, datetime
import csv


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
    # This file has no headers so here's an example instead:
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
    print(recent_pubmed_array[0])

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

def get_gene_info(species_gene_info_tsv, gene_citation_counts, taxonomy_name_tsv):
    """ 
    Getting a list of the gene information per gene and species' taxonomy name.
    """
    
    # Adding gene's information to the list created from get_recent_gene_citation_count per gene and saving it to the gene_info variable
    gene_symbol__gene_info__dict = {}
    # This file has no headers so here's an example instead:
    # 9606    59272   ACE2    -       ACEH    MIM:300335|HGNC:HGNC:13557|Ensembl:ENSG00000130234      X       Xp22.2  angiotensin converting enzyme 2 protein-coding  ACE2    angiotensin converting enzyme 2 O       angiotensin-converting enzyme 2|ACE-related carboxypeptidase|angiotensin I converting enzyme (peptidyl-dipeptidase A) 2|angiotensin I converting enzyme 2|angiotensin-converting enzyme homolog|angiotensin-converting enzyme-related carboxypeptidase|metalloprotease MPROT15|peptidyl-dipeptidase A|truncated angiotensin converting enzyme 2 20210711        -
    with open(species_gene_info_tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            if row[0][0] == '#':
                continue
            gene_id = str(row[1])
            tax_id = str(row[0])
            symbol = str(row[2]).upper()
            chromosome = str(row[6])
            full_name = str(row[8])
            gene_symbol__gene_info__dict[symbol] = {
                "symbol": symbol,
                "gene_id": gene_id,
                "chromosome": chromosome,
                "full_name": full_name,
                "citation_count": gene_citation_counts.get(gene_id, 0),
                "tax_id": tax_id
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
                        gene_symbol = extract_value_from_unparsed_string(details, "gene_id").upper()
                        start_coordinate = int(row[3])
                        end_coordinate = int(row[4])
                        coordinate_length = abs(end_coordinate - start_coordinate)
                        color = "#73af42"
                        if gene_symbol in gene_symbol__gene_info__dict.keys():
                            chromosome = gene_symbol__gene_info__dict[gene_symbol]["chromosome"]
                            full_name = gene_symbol__gene_info__dict[gene_symbol]["full_name"]
                            citation_count = gene_symbol__gene_info__dict[gene_symbol]["citation_count"]
                        else:
                            chromosome = ""
                            full_name = ""
                            citation_count = 0
                        ref_gene[gene_symbol] = {
                            "chromosome": chromosome,
                            "start_coordinate": start_coordinate,
                            "coordinate_length": coordinate_length,
                            "color": color,
                            "full_name": full_name,
                            "citation_count": citation_count,
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
    sorted_genes_list = sorted(ref_gene.items(), key=lambda x: x[1]['citation_count'], reverse=True)
    
    print("Gene with the most citations: " + str(sorted_genes_list[0]))

    return sorted_genes_list

def create_tsv_for_genes(sorted_genes_list, tax_name):
    """
    Creating TSVs containing gene information.
    """
    
    # Getting the species' scientific taxonomy name and saving it to the gene_species_name variable
    gene_species_name = tax_name.replace(" ", "-").replace("(", "-").replace(")", "-").replace("/", "-").replace("=", "-").lower()
    
    # Outputting the gene_species_name variable
    print(gene_species_name)

    # Getting today date
    today = date.today().strftime("%Y_%m_%d")

    # Setting the 5 months time frame
    days = 155

    # Getting the date 5 months ago
    month_ago = (datetime.today() - timedelta(days=155)).strftime("%Y_%m_%d")
    
    # Creating the citation variable that states the time frame from five months ago to today
    citations = f"citations_from_{month_ago}_to_{today}"
    
    # Creating the TSV name that contains the gene_species_name variable
    tsv_name= f'data/{gene_species_name}-citation-information.tsv'
    
    # Creating the TSV
    with open(tsv_name, 'wt') as out_file:
        # Setting up the TSV
        tsv_writer = csv.writer(out_file, delimiter='\t')
        
        # Add the header row to the TSV
        tsv_writer.writerow(["#name", "chromosome", "start", "length", "color", "full_name", citations, "significance"])
        
        # Going through each gene in the sorted_genes_list 
        for gene in sorted_genes_list:
            # Adding the gene information to the TSV
            symbol = gene[0]
            chromosome = gene[1]["chromosome"]
            start = gene[1]["start_coordinate"]
            length = gene[1]["coordinate_length"]
            color = gene[1]["color"]
            full_name = gene[1]["full_name"]
            citations = gene[1]["citation_count"]
            significance = "" # `significance` is now being handled from the front-end. leaving this here for backwards compatibility. Remove when no longer needed.
            tsv_row_values = [symbol, chromosome, start, length, color, full_name, citations, significance]

            tsv_writer.writerow(tsv_row_values)
    
# Main Function
if __name__ == "__main__":
    # Setting list of species
    list_of_species = ["human", "mouse", "rat"]

    # Going through each species
    for species in list_of_species:
        # Outputing species name
        print(species)
        
        gene_citation_counts = get_recent_gene_citation_count(f"creating_citation_counts_tsv/data/{species}/gene2pubmed", "creating_citation_counts_tsv/data/recent_pmid_year.ssv")
        
        gene_info, tax_name = get_gene_info(f"creating_citation_counts_tsv/data/{species}/gene_info", gene_citation_counts, "creating_citation_counts_tsv/taxonomy_name")
        
        ref_gene = get_ref_gene(species, gene_info)

        sorted_genes_list = sort_and_list_genes(ref_gene)

        create_tsv_for_genes(sorted_genes_list, tax_name)
