""" 
Creating TSVs containing the 10 most-cited genes per species and their gene information.
"""
import os
import sys
import json
import pyspark
import requests
from operator import add
from datetime import date, timedelta, datetime
import csv
from collections import Counter


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
    with open(species_gene_info_tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            if row[0][0] == '#':
                continue
            gene_id = str(row[1])
            tax_id = str(row[0])
            symbol = str(row[2]).upper()
            chromosome = str(row[6])
            description = str(row[8])
            gene_symbol__gene_info__dict[symbol] = {
                "symbol": symbol,
                "gene_id": gene_id,
                "chromosome": chromosome,
                "description": description,
                "citation_count": gene_citation_counts.get(gene_id, 0),
                "tax_id": tax_id
            }
    
    # Outputing the first element in the gene_info variable
    first_dict_entry = gene_symbol__gene_info__dict[list(gene_symbol__gene_info__dict.keys())[0]]
    print(first_dict_entry)
    
    # Getting the taxonomy ID for the species and saving it to the species_tax_id variable
    species_tax_id = first_dict_entry["tax_id"]

    # Outputing the species_tax_id variable
    print(species_tax_id)

    # Getting the taxonomy scientific name for the taxonomy ID given and saving it to the taxonomy_scientific_name variable
    taxonomy_scientific_name = ""
    # This file has no headers so here's an example instead:
    #28      |       halophilic eubacterium  |               |       scientific name |
    with open(taxonomy_name_tsv) as fd:
        rd = csv.reader(fd, delimiter='\t')
        for row in rd:
            row_tax_id = int(row[0])
            name = str(row[2]) # the file does not include column headers so these are guesses based on the data
            name_type = str(row[6])
            if row_tax_id==species_tax_id and name_type=="scientific name":
                taxonomy_scientific_name = name
                break

    # Outputing the taxonomy_scientific_name variable
    print(taxonomy_scientific_name)

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
                            description = gene_symbol__gene_info__dict[gene_symbol]["description"]
                            citation_count = gene_symbol__gene_info__dict[gene_symbol]["citation_count"]
                        else:
                            chromosome = ""
                            description = ""
                            citation_count = 0
                        ref_gene[gene_symbol] = {
                            "chromosome": chromosome,
                            "start_coordinate": start_coordinate,
                            "coordinate_length": coordinate_length,
                            "color": color,
                            "description": description,
                            "citation_count": citation_count,
                            }

    # Outputing the ref_gene variable
    print( next(iter( ref_gene.items() )) )
    
    # Returning the ref_gene variable
    return ref_gene

def translate_disease(disease):
    """ 
    Getting the common translated disease per disease given
    """
    
    # Dictionary of common translated disease to disease
    disease_translations = {
    'Neoplasms' : 'Cancer',
    'Musculoskeletal Diseases' : 'Muscle and Bone Diseases',
    'Digestive System Diseases' : 'Gastrointestinal Diseases',
    'Stomatognathic Diseases' : 'Mouth and Jaw Diseases',
    'Respiratory Tract Diseases' : 'Respiratory System Diseases',
    'Otorhinolaryngologic Diseases' : 'Ears, Throat and Nose Diseases',
    'Nervous System Diseases' : 'Nervous System Diseases',
    'Eye Diseases' : 'Eye Diseases',
    'Male Urogenital Diseases' : 'Male Genital Organ Diseases',
    'Female Urogenital Diseases and Pregnancy Complications' : 'Female Genital Organ and Pregnancy Diseases',
    'Cardiovascular Diseases' : 'Heart and Blood Diseases',
    'Hemic and Lymphatic Diseases' : 'Blood Diseases',
    'Congenital, Hereditary, and Neonatal Diseases and Abnormalities' : 'Inherited and Newborn Diseases',
    'Skin and Connective Tissue Diseases' : 'Skin Diseases',
    'Nutritional and Metabolic Diseases' : 'Nutrition and Metabolism Related Diseases',
    'Endocrine System Diseases' : 'Hormone Related Diseases',
    'Immune System Diseases' : 'Immune System Diseases',
    'Disorders of Environmental Origin' : 'Environmental Disorders',
    'Occupational Diseases' : 'Work-Related Diseases',
    'Substance-Related Disorders' : 'Substance-Related Disorders',
    'Wounds and Injuries' : 'Wounds and Injuries',
    'Behavior and Behavior Mechanisms' : 'Behavior',
  }
    try: 
        # Check to see if the disease is disease_translations Dictionary
        # Return common translated disease (If found)
        return disease_translations[disease]
    except:
        # Return 'Other' (If not found)
        return 'Other'

def request_disgenet(endpoint, params=None):
    """ 
    Makes a request to GWAS API with specified endpoint and parameters.
    Returns the response in json.
    Raises on request failure.
    """
    
    # Setting base_url variable to input in the Disgenet API
    base_url = "https://www.disgenet.org/api"
    
    # Setting base_url variable to input in the Disgenet API
    headers = {"content-type": "application/json"}
    
    # Calling the Disgenet API
    try:
        # Calling the Disgenet API and saving it in the response variable
        response = requests.get(
            url=f"{base_url}/{endpoint}", params=params, headers=headers
        )
        # Check the response variable status
        response.raise_for_status()
        # Return the response variable as a JSON
        return response.json()
    except (
        # Rasing errors the Disgenet API call fails
        requests.exceptions.ConnectionError,
        requests.exceptions.HTTPError,
        requests.exceptions.Timeout,
    ) as e:
        # Raise RuntimeError(f"Request failed: {e}")
        return [{"disease_class_name":[]}]

def get_significance(gene, source, num_diseases):
    """
    Getting all diseases associated to gene.
    """
    # Create gda (gene-disease association) variable to input in the Disgenet API
    gda = "gda/gene/"
    
    # API call to Disgenet database to get a list of diseases for the gene
    gene_related_diseases = request_disgenet(endpoint=f"{gda}{gene}", params=f"source={source}&format=json")
    
    # Creating the disease_associations variable
    disease_associations = []
    
    # Going through each disease in the gene_related_diseases list 
    for disease in gene_related_diseases:
        # Skipping instances with no disease class name
        if disease['disease_class_name']:
            # Getting base disease class and strip whitespaces
            disease_associations.append(disease['disease_class_name'].split(';')[0].strip())
  
    # Checks if disease_associations list is empty 
    if not disease_associations:
        # Adding None to the disease_associations List
        disease_associations.append("None")

    # Talling disease class names
    tallied_disease_associations = Counter(disease_associations)
  
    # Getting top n + 3 gene disease associations
    # + 3 to account for multiple 'other' disease classes (update to n-related variable)
    disease_associations = [association for (association, _) in tallied_disease_associations.most_common(num_diseases + 3)]
    
    # Creating the top_disease_associations variable
    top_disease_associations = []
    
    # Going through each disease in the disease_associations list 
    for disease in disease_associations: 
        # Limiting top disease associations
        if len(top_disease_associations) < num_diseases: 
            # Overwriting with layperson understandable terms 
            disease = translate_disease(disease) 
        # Allowing single other association
        if disease != 'Other' or 'Other' not in top_disease_associations:
            # Adding disease to top_disease_associations
            top_disease_associations.append(disease)

    # Check if other is in the top_disease_associations list
    if 'Other' in top_disease_associations:
        # Move Other to the back of the list
        top_disease_associations.append(top_disease_associations.pop(top_disease_associations.index('Other')))

    # Format and return top_disease_associations
    return f"Involved in {'; '.join(top_disease_associations).lower()}" if top_disease_associations else None


def create_tsv_for_most_cited_genes(ref_gene, tax_name, significance_SOURCES, top_count):
    """
    Creating TSVs containing the `top_count` most-cited genes per species and their gene information.
    """
    
    print("Getting Top " + str(top_count) + " most-cited Genes")
    
    # Getting the list of the top ten most cited genes from the list created from get_ref_gene
    top_genes_list = sorted(ref_gene.items(), key=lambda x: x[1]['citation_count'], reverse=True)[:top_count]
    
    print("Gene with the most citations:")
    print(str(top_genes_list[0]))
    
    # Creating the top_genes_with_significance_list variable
    top_genes_with_significance_list = []
    
    # Going through each gene in the top_genes_list 
    for gene_row in top_genes_list:
                
        # Getting the top three name disease and saving it to the significance variable
        gene_name = gene_row[0]
        significance = get_significance(gene_name, significance_SOURCES, 3)
        
        # Add the significance information to the gene_row list
        gene_row = gene_row + tuple([significance])
        
        # Adding the updated gene_row to the top_genes_with_significance_list
        top_genes_with_significance_list.append(gene_row)
    
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
    tsv_name= f'tsv/{gene_species_name}-citation-information.tsv'
    
    # Creating the TSV
    with open(tsv_name, 'wt') as out_file:
        # Setting up the TSV
        tsv_writer = csv.writer(out_file, delimiter='\t')
        
        # Add the header row to the TSV
        tsv_writer.writerow(["#name", "chromosome", "start", "length", "color", "full_name", citations, "significance"])
        
        # Going through each gene in the top_genes_with_significance_list 
        for gene in top_genes_with_significance_list:
            # Adding the gene information to the TSV
            tsv_writer.writerow(gene)
    
# Main Function
if __name__ == "__main__":
    # Setting list of species
    list_of_species = ["human", "mouse", "rat"]
    
    # Setting dictionary of species to Disgenet API source 
    significance_SOURCES = {
        "human": 'CTD_human' , 
        "mouse": 'MGD',
        "rat": 'RGD',
    }
    
    # Going through each species
    for species in list_of_species:
        # Outputing species name
        print(species)
        
        gene_citation_counts = get_recent_gene_citation_count(f"creating_citation_counts_tsv/data/{species}/gene2pubmed", "creating_citation_counts_tsv/data/recent_pmid_year.ssv")
        
        gene_info, tax_name = get_gene_info(f"creating_citation_counts_tsv/data/{species}/gene_info", gene_citation_counts, "creating_citation_counts_tsv/taxonomy_name")
        
        ref_gene = get_ref_gene(species, gene_info)

        create_tsv_for_most_cited_genes(ref_gene, tax_name, significance_SOURCES[species], 10)
        #todo: remove citations counts == 0

        
#todo -- add csv comments to each file open
