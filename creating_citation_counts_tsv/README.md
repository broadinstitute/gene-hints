# Creating Citation Count TSV 

![TSV_snapshot](https://github.com/broadinstitute/gene-hints/blob/getting_tsv/images/07-TSV-snapshot.png)

## Background: 
The Gene Hint UI uses data from generated TSVs to calculate the current top genes per species. This code will use the Entrez database to collect all the gene's publication citations over the course of five months. Next, the code will use NCBI's database to match the gene IDs to their gene publication. The code will then use the UCSC and GenBank-data database to match the gene IDs to all the gene's information (ex. Chromosome and position). The code will then use the DisGeNet database to match the gene ID to the gene's disease associations (ex. cancer or heart disease). Lastly, the code will output TSVs, of the 10 most cited genes and their gene information.

## To Run: 
1. In the root directory of gene-hints, run the create_citations_tsv shell script:

    ``` 
    sh creating_citation_counts_tsv/scripts/create_citations_tsv.sh 
    ```

## Adding More Species:

Currently, this code will produce TSVs human, mouse, and rat species. To add another species, you will need to add 4 lines of code to the create_citations_tsv.sh shell script and make two adjustments to the summarize_gene_citations_all_species.py python scripts. 

1. ###  Additions to create_citations_tsv.sh  
    1.1. Add a line that will create new directory for your species **(Line 35)**
    
     ``` 
     mkdir creating_citation_counts_tsv/data/{species_name} 
     ```
        
     1.2. Each species has one or more taxonomy IDs. You can look on this website https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi to find the tax ID for your species. If you want all the tax IDs regarding your species, which we suggest, click the full tax ID lineage check box before searching. Once you find the tax ID, add a line that will sort through the genes to the publication ID file for your species and save it in your species' directory **(Line 55)**
    
      ``` 
       cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == {tax_ID_1} || $1 == {tax_ID_2) print;}' > creating_citation_counts_tsv/data/{species_name}/gene2pubmed
      ```
      1.3. Each species have a specific reference file. Using this website, https://hgdownload.soe.ucsc.edu/downloads.html, you need to find the URL for your species' reference gene file ending in ".gtf.gz". After you found the refGene URL, add a line that will download your species' reference gene.**(Line 71)**
    
      ``` 
        wget -N -P creating_citation_counts_tsv/data/{species_name} 
        
        gunzip creating_citation_counts_tsv/data/{species_name}/{file_name}
      ```
        
      1.4. Using the same tax IDs, add a line that will sort through the gene info file for your species and save it in your species' directory **(Line 84)**
    
      ``` 
        cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == {tax_ID_1} || $1 == {tax_ID_2}) print;}' > creating_citation_counts_tsv/data/{species_name}/gene_info 
      ```
2. ###  Additions to summarize_gene_citations_all_species.py
    2.1. In the list_of_species list in the main function, add your species to the list
    
     ``` 
        list_of_species = ["human", "mouse", "rat", {species_name} ]
     ```
     2.2. For the DisGeNet database to find the related diseases per gene, you have to pick a model source, so your genes will be associated with the right species. This (https://www.disgenet.org/api/#/GDA/gdaBySource) will display all the sources you can choose from. Pick which source relates to your species. Then you can add a line that will add your species name and your species source to the significance_SOURCES dictionary.
    
     ``` 
        significance_SOURCES = {
        "human": 'CTD_human' , 
        "mouse": 'MGD',
        "rat": 'RGD',
        "{species_name}": '{species_source}',
        }
     ```

