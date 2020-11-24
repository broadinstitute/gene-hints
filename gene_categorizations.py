'''
Get gene disease associations using the disgenet API
https://www.disgenet.org/api/
'''

import requests
import json
import sys

from collections import Counter

base_url = "https://www.disgenet.org/api"

# gene-disease associations
gda = "gda/gene/"

pubmed_id_search = "search/findByPubmedId?pubmedId="

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
        raise RuntimeError(f"Failed to connect to GWAS API: {e}")

def get_diseases():
    make_request(endpoint=f"{gda}")

def get_gene_disease_associations(genes):
    # wrap single gene in a list so we can handle single and list input types
    if not isinstance(genes, list):
        genes = [genes]

    # compile dictionary of pubmed IDs to EFO traits
    GDA = {}
    for gene in genes: 
        gene_related_diseases = make_request(endpoint=f"{gda}{gene}")
        GDA[gene] = Counter([disease['disease_semantic_type'] for disease in gene_related_diseases])
        
    print(GDA)
    return GDA
   
if len(sys.argv) < 1:
    sys.exit("USAGE: "
             "python3 gene_categorizations.py pubmedID"
             "where pubmedID can be a single ID or a list of IDs")

genes = sys.argv[1]

get_gene_disease_associations(genes)

'''
HTTP/1.1 200 OK
Content-Type: application/json;charset=UTF-8
Content-Length: 747

{
  "_links" : {
    "efoTraits" : {
      "href" : "https://www.ebi.ac.uk/gwas/rest/api/efoTraits{?page,size,sort,projection}",
      "templated" : true
    },
    "associations" : {
      "href" : "https://www.ebi.ac.uk/gwas/rest/api/associations{?page,size,sort,projection}",
      "templated" : true
    },
    "singleNucleotidePolymorphisms" : {
      "href" : "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms{?page,size,sort,projection}",
      "templated" : true
    },
    "studies" : {
      "href" : "https://www.ebi.ac.uk/gwas/rest/api/studies{?page,size,sort,projection}",
      "templated" : true
    },
    "profile" : {
      "href" : "https://www.ebi.ac.uk/gwas/rest/api/profile"
    }
  }
}
'''