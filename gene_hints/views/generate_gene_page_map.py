from SPARQLWrapper import SPARQLWrapper, JSON
from pandas import json_normalize
import os

output_location = "./wikipedia_trends/"
filename = "gene_page_map.tsv"


def query_wikidata(sparql_query, sparql_service_url):
    """Query endpoint with given query string and return the results as a
    pandas Dataframe.
    """
    sparql = SPARQLWrapper(sparql_service_url, agent="chrome")

    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)

    result = sparql.query().convert()
    return json_normalize(result["results"]["bindings"])


def query_human_genes():
    """Execute the SPARQL query for human genes
    """
    print("Querying human genes...")
    endpoint_url = "https://query.wikidata.org/sparql"
    query = """
    SELECT DISTINCT ?item ?ncbi_gene ?itemLabel ?titleLabel WHERE {
        ?gene schema:about ?item;
            schema:isPartOf <https://en.wikipedia.org/>;
            schema:name ?title.
        BIND(REPLACE(STR(?title), "\\\\ ", "_") AS ?titleLabel)
            ?item wdt:P351 ?ncbi_gene;
            wdt:P703 wd:Q15978631.
        SERVICE wikibase:label
        { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
    }"""

    return query_wikidata(query, endpoint_url)


def save_human_genes(data, name):
    """Save results of the gene locally
    """
    print("Saving results of gene locally...")
    data[["titleLabel.value", "itemLabel.value"]].rename(
        columns=lambda col: col.replace("Label.value", "")
    ).to_csv(name, sep="\t", index=False)


genes = query_human_genes()
save_human_genes(genes, os.path.join(output_location, filename))
