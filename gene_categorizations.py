'''
Get gene disease associations using the disgenet API
https://www.disgenet.org/api/
'''

import requests
import json
import csv
import sys
import os 

from collections import Counter

base_url = "https://www.disgenet.org/api"

# gene-disease associations
gda = "gda/gene/"
disease = "disease/"

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
  
  # get top n gene disease associations
  top_disease_associations = [association for (association, _) in tallied_disease_associations.most_common(num_diseases)]

  # overwrite with layperson understandable terms
  top_disease_associations = [translate_disease(disease) for disease in top_disease_associations]
  
  # format and return
  return f"Involved in {'; '.join(top_disease_associations).lower()}"

SOURCES = {
  'homo-sapiens': 'CTD_human' , 
  'mus-musculus': 'MGD',
  'rattus-norvegicus': 'RGD',
}

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
              gene_name = row[0].upper()
              
              if '#NAME' in gene_name:
                # copy over the header 
                processed_row = row
              else:
                try:
                  # get each gene's significance
                  significance = get_significance(gene_name, SOURCES[source], num_diseases)
                  row[7] = significance
                  processed_row = row
                except:
                  print(f"Failed to retrieve significance, skipping {gene_name}.")
                  continue
              # write to the processed file
              writer.writerow(processed_row) 
                
# run the script! 
add_significance_to_TSVs()

