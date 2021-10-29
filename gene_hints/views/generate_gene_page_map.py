"""Output TSV of English Wikipedia article names to gene symbols

For example, the page https://en.wikipedia.org/wiki/Tumor_necrosis_factor maps
to the gene symbol "TNF".  This is output as:

Tumor_necrosis_factor	TNF

in `gene_page_map.tsv`.  This script uses the Wikidata Query Service
(https://query.wikidata.org/) to make a SPARQL query linking articles to
genes symbols.  The output TSV is used by `views.py`.
"""

from SPARQLWrapper import SPARQLWrapper, JSON
from pandas import json_normalize

output_path = "./data/gene_page_map.tsv"

def query_wikidata(sparql_query, sparql_service_url):
    """Query endpoint with given query string and return the results as a
    pandas Dataframe.
    """
    sparql = SPARQLWrapper(sparql_service_url, agent="chrome")

    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)

    result = sparql.query().convert()
    return json_normalize(result["results"]["bindings"])

def query_human_genes_special():
    """Execute SPARQL query for article names for special human genes

    Some articles about genes are modeled as proteins on Wikidata.
    These tend to be prominent genes, like BRCA1 and TP53.
    """
    print("Querying human genes, special case...")
    endpoint_url = "https://query.wikidata.org/sparql"
    query = """
    SELECT DISTINCT ?item ?gene_symbol ?titleLabel WHERE {
        ?protein schema:about ?item;
                schema:isPartOf <https://en.wikipedia.org/>;
                schema:name ?title.
        BIND(REPLACE(STR(?title), "\\\\ ", "_") AS ?titleLabel)
            ?item wdt:P31 wd:Q8054. # is a protein
            ?item wdt:P703 wd:Q15978631. # found in human
            ?item wdt:P702 ?gene.
            ?gene wdt:P353 ?gene_symbol.
    }
    """
    return query_wikidata(query, endpoint_url)

def query_human_genes_general():
    """Execute SPARQL query for article names for general human genes

    Almost all articles about genes are modeled as genes on Wikidata.
    This handles the common case.
    """
    print("Querying human genes, general case...")
    endpoint_url = "https://query.wikidata.org/sparql"
    query = """
    SELECT DISTINCT ?item ?ncbi_gene ?itemLabel ?titleLabel WHERE {
        ?gene schema:about ?item;
            schema:isPartOf <https://en.wikipedia.org/>;
            schema:name ?title.
        BIND(REPLACE(STR(?title), "\\\\ ", "_") AS ?titleLabel)
            ?item wdt:P351 ?ncbi_gene;
            wdt:P703 wd:Q15978631. # found in human
        SERVICE wikibase:label
        { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
    }"""

    return query_wikidata(query, endpoint_url)

def save_human_genes(data_special, data_general, output_path):
    """Save results of the gene query locally
    """
    print("Saving results of gene query locally...")
    data_special[["titleLabel.value", "gene_symbol.value"]].rename(
        columns=lambda col: col.replace("Label.value", "")
    ).to_csv(output_path, sep="\t", index=False)
    data_general[["titleLabel.value", "itemLabel.value"]].rename(
        columns=lambda col: col.replace("Label.value", "")
    ).to_csv(output_path, sep="\t", index=False, mode="a", header=False)
    print("Results saved to: " + output_path)


genes_special = query_human_genes_special()
genes_general = query_human_genes_general()
save_human_genes(genes_special, genes_general, output_path)
