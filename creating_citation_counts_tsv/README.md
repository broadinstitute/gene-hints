# Creating Citation Count TSV 

![TSV_snapshot](https://github.com/broadinstitute/gene-hints/blob/getting_tsv/images/07-TSV-snapshot.png)

## Background: 
The Gene Hint UI uses data from generated TSVs to calculate the current top genes per species. This code will use the Entrez database to collect all the gene's publication citations over the course of five months. Next, the code will use NCBI's database to match the gene IDs to their gene publication. The code will then use the UCSC and GenBank-data database to match the gene IDs to all the gene's information (ex. Chromosome and position). The code will then use the Disgenet database to match the gene ID to the gene's disease associations (ex. cancer or heart disease). Lastly, the code will output TSVs, of the 10 most cited genes and their gene information.

## Setup: 

## To Run: 
1. In the root directory of gene-hints, run the create_citations_tsv shell script:

    ``` 
    sh creating_citation_counts_tsv/scripts/create_citations_tsv.sh 
    ```

## Adding More Species: 
