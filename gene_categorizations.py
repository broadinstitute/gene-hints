'''
Get gene disease associations using the disgenet API
https://www.disgenet.org/api/
'''

import requests
import json
import csv
import sys
import os 

from itertools import takewhile
from collections import Counter

base_url = "https://www.disgenet.org/api"

# gene-disease associations
gda = "gda/gene/"
disease = "disease/"

SOURCES = {
  'homo-sapiens': 'CTD_human' , 
  'mus-musculus': 'MGD',
  'rattus-norvegicus': 'RGD',
}


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
        raise RuntimeError(f"Request failed: {e}")

def get_significance(gene, source, num_diseases): 
  # get all diseases associated to gene
  gene_related_diseases = make_request(endpoint=f"{gda}{gene}", params=f"source={source}&format=json")
  disease_associations = []
  for disease in gene_related_diseases:
    # skip instances with no disease class name
    if disease['disease_class_name']:
      # get base disease class and strip whitespaces
      disease_associations.append(disease['disease_class_name'].split(';')[0].strip())
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
  return f"Involved in {'; '.join(top_disease_associations).lower()}"


def add_significance_to_TSVs():
  num_diseases = 3
  original_dir = 'all_species_tsv'
  processed_dir = 'processed_tsv'

  for filename in os.listdir(original_dir):
    # only process human, rat and mouse TSV files
    for source in SOURCES.keys():
      if filename.startswith(source) and filename.endswith(".tsv"): 
        print(f"Processing {source} file.")
        # set up file reader
        with open(f"{original_dir}/{filename}", encoding = "ISO-8859-1") as read_file:
          reader = csv.reader(read_file, delimiter="\t")
          
          # set up file writer
          with open(f"{processed_dir}/{filename}", 'wt') as write_file:
            writer = csv.writer(write_file, delimiter='\t')
            for row in reader:
              gene_name = row[0]
              
              # skip significance search for header
              if '#name' not in gene_name:
                try:
                  # get each gene's significance
                  significance = get_significance(gene_name, SOURCES[source], num_diseases)
                  row[7] = significance
                except:
                  print(f"Failed to retrieve significance, skipping {gene_name}.")
                
              # write to the processed file
              processed_row = row
              writer.writerow(processed_row) 
                
# run the script! 
add_significance_to_TSVs()

