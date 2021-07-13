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


def get_recent_gene2pubmed(species):
    """ 
    Getting a list of the gene ID, gene publication citation ID, and count of how frequently the gene is cited over the past 5 months for each gene.
    """

    # Create a dict mapping pubmed_id to gene_id
    pubmed_gene_dict = {}
    with open(f"creating_citation_counts_tsv/data/{species}/gene2pubmed") as fd:
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
    with open("creating_citation_counts_tsv/data/recent_pmid_year.ssv") as fd:
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
    gene_counts = {}
    for gene_id in recent_gene_pubmed_dict:
        gene_counts[gene_id] = len(recent_gene_pubmed_dict[gene_id])

    # Outputing the first element in the gene_counts variable
    print( next(iter( gene_counts.items() )) )
    
    # Returning the gene_counts variable
    return gene_counts

def getting_gene_info(species, gene_counts):
    """ 
    Getting a list of the gene information per gene and species' taxonomy name.
    """
    
    # Adding gene's information to the list created from get_recent_gene2pubmed per gene and saving it to the gene_info variable
    gene_info = (sc.textFile(f"creating_citation_counts_tsv/data/{species}/gene_info")
                    .filter(lambda x: x[0] != '#')
                    .map(lambda x: x.split('\t'))
                    .map(lambda x: (str(x[1]), (str(x[0]), str(x[2]), str(x[6]), str(x[8]))))).join(gene_counts).map(lambda x: (x[1][0][1].upper(), (x[0], x[1][0][2], x[1][0][3], x[1][1][1], x[1][0][0])))
    
    # Outputing the first element in the gene_info variable
    print(gene_info.take(1))
    
    # Getting the taxonomy ID for the species and saving it to the tax_id variable
    tax_id = str(gene_info.collect()[0][1][4])
    
    # Outputing the first element in the tax_id variable
    print(tax_id)

    # Getting the taxonomy information for the taxonomy ID given and saving it to the taxonomy variable
    taxonomy = (sc.textFile("creating_citation_counts_tsv/taxonomy_name")
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .filter(lambda x: int(tax_id) == int(x[0]) )
                        .filter(lambda x: "scientific name" in x[6] )
                        .map(lambda x: ((str(x[0]), str(x[2]), str(x[6]))))
                     )
    
    # Getting the taxonomy name from the taxonomy information and saving it to the tax_name variable
    tax_name = taxonomy.collect()[0][1]
    
    # Outputing the tax_name variable
    print(tax_name)

    # Returning the gene_info and tax_name variable
    return gene_info, tax_name

def get_ref_gene(species, gene_info): 
    """ 
    Getting a list of the gene reference information per gene.
    """
    
    # Going through all the files in the species directory 
    for refGene_file in os.listdir(f"creating_citation_counts_tsv/data/{species}"):
         #Since the reference file is uniquely named and could not be hard-coded, this will filter out "gene_info" and "gene2pubmed" to read the reference file
        if refGene_file != "gene_info" and "gene2pubmed":
            # Adding gene's reference information to the list created from getting_gene_info per gene and saving it to the refGene variable
            refGene = (sc.textFile(f"creating_citation_counts_tsv/data/{species}/{refGene_file}")
                    .filter(lambda x: x[0] != '#')
                    .map(lambda x: x.split('\t'))
                    .filter(lambda x: 'transcript' in x[2])
                    .map(lambda x: (str(x[8]).split(";")[0].replace("gene_id ", "").replace('"', '').upper(), (int(x[3]), int(abs(int(x[3])-int(x[4])))))).join(gene_info).map(lambda x: (x[0], (x[1][1][1], x[1][0][0], x[1][0][1], "#73af42", x[1][1][2], x[1][1][3]))).reduceByKey(lambda a, b: a).map(lambda x: (x[0], x[1][0], x[1][1], x[1][2], "#73af42", x[1][4], x[1][5])))
            
            # Outputing the refGene variable
            print(refGene.take(1))
            
            # Returning the refGene variable
            return refGene

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


def create_tsv(refGene, tax_name, significance_SOURCES):
    """
    Creating TSVs containing the 10 most-cited genes per species and their gene information.
    """
    
    # Outputing the start of the Getting Top Ten TSVs process
    print("Getting Top Ten TSVs")
    
    # Getting the list of the top ten most cited genes from the list created from get_ref_gene
    top_ten_gene = sorted(refGene.collect(), key=lambda x: x[6])[-10:][::-1]
    
    # Creating the top_ten_gene_list variable
    top_ten_gene_list = []
    
    # Going through each gene in the top_ten_gene list 
    for gene_row in top_ten_gene:
        # Getting the gene name and saving it to the gene_name variable
        gene_name = gene_row[0]
        
        # Getting the top three name disease and saving it to the significance variable
        significance = get_significance(gene_name, significance_SOURCES, 3)
        
        # Add the significance information to the gene_row list
        gene_row = gene_row + tuple([significance])
        
        # Adding the updated gene_row to the top_ten_gene_list
        top_ten_gene_list.append(gene_row)
    
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
        
        # Going through each gene in the top_ten_gene_list 
        for gene in top_ten_gene_list:
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
        
        # Getting a list of the gene ID, gene publication citation ID, and count of how frequently the gene is cited over the past 5 months for each gene
        gene_counts = get_recent_gene2pubmed(species)
        
        # Getting a list of the gene information per gene and species' taxonomy name.
        gene_info, tax_name = getting_gene_info(species, gene_counts)
        
        # Getting a list of the gene reference information per gene
        refGene = get_ref_gene(species, gene_info)
        
        # Creating TSVs containing the 10 most-cited genes per species and their gene information
        create_tsv(refGene, tax_name, significance_SOURCES[species])

        
