'''
Get gene associations for a pubmed ID or list of pubmed IDs
https://www.ebi.ac.uk/gwas/rest/docs/api#resources-associations
'''

import requests
import json
import sys

base_url = "https://www.ebi.ac.uk/gwas/rest/api/"

studies = "studies/"
associations = "associations/"

pubmed_id_search = "search/findByPubmedId?pubmedId="
study_associations = "associations?projection=associationsByStudySummary"

def make_gwas_request(endpoint, params=None):
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

def get_associations(pubmed_ids):
    # wrap single pubmed IDs in a list so we can handle single and list input types
    if not isinstance(pubmed_ids, list):
        pubmed_ids = [pubmed_ids]

    # compile dictionary of pubmed IDs to associations
    PM_associations = {}
    for pubmed_id in pubmed_ids: 
        data = make_gwas_request(endpoint=studies, params=f"{pubmed_id_search}{pubmed_id}{study_associations}")["_embedded"]["studies"]
        print(data)
        #PM_associations = data['association']
        #PM_associations[pubmed_id] = association
    
    #return PM_associations

if len(sys.argv) < 1:
    sys.exit("USAGE: "
             "python3 gene_categorizations.py pubmedID"
             "where pubmedID can be a single ID or a list of IDs")

pubmed_ids = sys.argv[1]

get_associations(pubmed_ids)

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