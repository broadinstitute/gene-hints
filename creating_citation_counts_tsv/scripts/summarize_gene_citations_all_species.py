import os
import findspark
import pyspark
import shortuuid
import requests
from tqdm import tqdm
from operator import add
from datetime import date, timedelta, datetime
import shutil
import time
import csv
from itertools import takewhile
from collections import Counter
import json
import sys

def getting_recent_gene2pubmed(species):
    gene2pubmed = (sc.textFile(f"creating_citation_counts_tsv/data/{species}/gene2pubmed")
                    .filter(lambda x: x[0] != '#')
                    .map(lambda x: x.split('\t'))
                    .map(lambda x: ( str(x[2]), int(x[1]))))
    print(gene2pubmed.take(1))
    t1 = time.time()
    pubmeds_set = set([x[0] for x in gene2pubmed.collect()])
    t2 = time.time()

    print("time taken", t2 - t1)

    year_pmid = (sc.textFile("creating_citation_counts_tsv/data/recent_pmid_year.ssv")
                    .map(lambda x: x.split())
                    .map(lambda x: (str(x[1]), int(x[0])))
                    .filter(lambda x: x[0] in pubmeds_set))
    print(year_pmid.take(1))
    
    gene2pubmed_to_year = gene2pubmed.join(year_pmid).map(lambda x: (x[1][0], x[0]))
    print(gene2pubmed_to_year.take(1))
    gene_counts = gene2pubmed.join(year_pmid).map(lambda x: (x[1][0], 1)).reduceByKey(add).join(gene2pubmed_to_year).map(lambda x: ( str(x[0]), (x[1][1] ,x[1][0])))
    print(gene_counts.take(1))
    return gene_counts

def getting_gene_info(species, gene_counts):
    gene_info = (sc.textFile(f"creating_citation_counts_tsv/data/{species}/gene_info")
                    .filter(lambda x: x[0] != '#')
                    .map(lambda x: x.split('\t'))
                    .map(lambda x: (str(x[1]), (str(x[0]), str(x[2]), str(x[6]), str(x[8]))))).join(gene_counts).map(lambda x: (x[1][0][1].upper(), (x[0], x[1][0][2], x[1][0][3], x[1][1][1], x[1][0][0])))
    print(gene_info.take(1))
    print(gene_info.collect()[0][1][4])
    tax_id = str(gene_info.collect()[0][1][4])
    print(tax_id)

    taxonomy = (sc.textFile("creating_citation_counts_tsv/taxonomy_name")
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .filter(lambda x: int(tax_id) == int(x[0]) )
                        .filter(lambda x: "scientific name" in x[6] )
                        .map(lambda x: ((str(x[0]), str(x[2]), str(x[6]))))
                     )
    print(taxonomy.collect()[0][1])
    tax_name = taxonomy.collect()[0][1]

    return gene_info, tax_name

def getting_refGene(species, gene_info):    
    print(f"Getting gene refGene for {species}")  
    for refGene_file in os.listdir(f"creating_citation_counts_tsv/data/{species}"):
        if refGene_file != "gene_info" and "gene2pubmed":
            refGene = (sc.textFile(f"creating_citation_counts_tsv/data/{species}/{refGene_file}")
                    .filter(lambda x: x[0] != '#')
                    .map(lambda x: x.split('\t'))
                    .filter(lambda x: 'transcript' in x[2])
                    .map(lambda x: (str(x[8]).split(";")[0].replace("gene_id ", "").replace('"', '').upper(), (int(x[3]), int(abs(int(x[3])-int(x[4])))))).join(gene_info).map(lambda x: (x[0], (x[1][1][1], x[1][0][0], x[1][0][1], "#73af42", x[1][1][2], x[1][1][3]))).reduceByKey(lambda a, b: a).map(lambda x: (x[0], x[1][0], x[1][1], x[1][2], "#73af42", x[1][4], x[1][5])))
            print(refGene.take(1))
            return refGene

def translate_disease(disease):
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
    return disease_translations[disease]
  except:
    return 'Other'

def make_request(endpoint, params=None):
    """ 
    Makes a request to GWAS API with specified endpoint and parameters.
    Returns the response in json.
    Raises on request failure.
    """
    base_url = "https://www.disgenet.org/api"
    headers = {"content-type": "application/json"}
    try:
        response = requests.get(
            url=f"{base_url}/{endpoint}", params=params, headers=headers
        )
        response.raise_for_status()
        return response.json()
    except (
        requests.exceptions.ConnectionError,
        requests.exceptions.HTTPError,
        requests.exceptions.Timeout,
    ) as e:
        # raise RuntimeError(f"Request failed: {e}")
        return [{"disease_class_name":[]}]

def get_significance(gene, source, num_diseases): 
  # get all diseases associated to gene
  gda = "gda/gene/"
  gene_related_diseases = make_request(endpoint=f"{gda}{gene}", params=f"source={source}&format=json")
  disease_associations = []
  for disease in gene_related_diseases:
    # skip instances with no disease class name
    if disease['disease_class_name']:
      # get base disease class and strip whitespaces
      disease_associations.append(disease['disease_class_name'].split(';')[0].strip())
  
  if not disease_associations:
    disease_associations.append("None")

  # tally disease class names
  tallied_disease_associations = Counter(disease_associations)
  
  # get top n + 3 gene disease associations
  # + 3 to account for multiple 'other' disease classes (update to n-related variable)
  disease_associations = [association for (association, _) in tallied_disease_associations.most_common(num_diseases + 3)]

  top_disease_associations = []
  for disease in disease_associations: 
    # limit top disease associations
    if len(top_disease_associations) < num_diseases: 
      # overwrite with layperson understandable terms 
      disease = translate_disease(disease) 
      # only allow single other association
      if disease != 'Other' or 'Other' not in top_disease_associations:
        top_disease_associations.append(disease)

  if 'Other' in top_disease_associations:
    # move Other to the back of the list
    top_disease_associations.append(top_disease_associations.pop(top_disease_associations.index('Other')))

  # format and return
  return f"Involved in {'; '.join(top_disease_associations).lower()}" if top_disease_associations else None


def create_tsv(refGene, tax_name, significance_SOURCES):
    print("Get Top Ten")
    top_ten_gene = sorted(refGene.collect(), key=lambda x: x[6])[-10:][::-1]
    top_ten_gene_list = []
    for gene_row in top_ten_gene:
        gene_name = gene_row[0]
        significance = get_significance(gene_name, significance_SOURCES, 3)
        gene_row = gene_row + tuple([significance])
        top_ten_gene_list.append(gene_row)
    
    gene_species_name = tax_name.replace(" ", "-").replace("(", "-").replace(")", "-").replace("/", "-").replace("=", "-").lower()
    print(gene_species_name)

    
    # Getting today date
    today = date.today().strftime("%Y_%m_%d")

    # Setting the 5 months time frame
    days = 155

    # Getting the date 5 months ago
    month_ago = (datetime.today() - timedelta(days=155)).strftime("%Y_%m_%d")
    citations = f"citations_from_{month_ago}_to_{today}"
    tsv_name= f'tsv/{gene_species_name}-citation-information.tsv'
    with open(tsv_name, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(["#name", "chromosome", "start", "length", "color", "full_name", citations, "significance"])
        for gene in top_ten_gene_list:
            tsv_writer.writerow(gene)
    

if __name__ == "__main__":
    list_of_species = ["human", "mouse", "rat"]
    significance_SOURCES = {
        "human": 'CTD_human' , 
        "mouse": 'MGD',
        "rat": 'RGD',
    }
    sc = pyspark.SparkContext()
    for species in list_of_species:
        print(species)
        gene_counts = getting_recent_gene2pubmed(species)
        gene_info, tax_name = getting_gene_info(species, gene_counts)
        refGene = getting_refGene(species, gene_info)
        create_tsv(refGene, tax_name, significance_SOURCES[species])

        