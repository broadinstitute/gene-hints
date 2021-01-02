#!/bin/bash

# Creating Citiation TSVs Per selected species. Adapted from source code: https://github.com/pkerpedjiev/gene-citation-counts

# First install 
pip3 install -r creating_citation_counts_tsv/requirements.txt 

# Create a data, tsv, and tmp-ssv folder
mkdir creating_citation_counts_tsv/data
mkdir tsv
mkdir creating_citation_counts_tsv/data/tmp-ssv

#To get an idea of how citation counts have changed over the years, 
#we need a list of which citations were published when. 
# This script queries Entrez gene for citation lists for each day for the time frame of 5 months.
python3 creating_citation_counts_tsv/scripts/pmids_by_date.py

# Consolidate the publications from each day into one complete list
# For easier processing, we'll collapse the per-day list of files into one file
for file in $(find . -name "*.ssv"); \
do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> creating_citation_counts_tsv/data/recent_pmid_year.ssv; \
done

# Remove tmp-ssv to save space
rm creating_citation_counts_tsv/data/tmp-ssv/*
rmdir creating_citation_counts_tsv/data/tmp-ssv

# Creating taxonomy folder 
mkdir creating_citation_counts_tsv/taxonomy

# Creating Species folders (human, mouse, rat)
mkdir creating_citation_counts_tsv/data/human
mkdir creating_citation_counts_tsv/data/mouse
mkdir creating_citation_counts_tsv/data/rat
# ADD Extra Species folders

# Getting taxonomy names
wget -N -P creating_citation_counts_tsv/taxonomy/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z 
zcat creating_citation_counts_tsv/taxonomy/new_taxdump.tar.Z > creating_citation_counts_tsv/taxonomy/new_taxdump.tar; /
tar -C creating_citation_counts_tsv/taxonomy/ -xvf creating_citation_counts_tsv/taxonomy/new_taxdump.tar
mv creating_citation_counts_tsv/taxonomy/names.dmp creating_citation_counts_tsv/taxonomy_name
rm creating_citation_counts_tsv/taxonomy/*

# The NCBI maintains a list of citations that mention each gene in its database. We'll download it.
wget -N -P creating_citation_counts_tsv/data ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz
gunzip creating_citation_counts_tsv/data/gene2pubmed.gz

# gene2pubmed is a big file, so we are parsing out what we need
# Getting Human gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 9606 || $1 == 9605) print;}' > creating_citation_counts_tsv/data/human/gene2pubmed
# Getting Mouse gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if (10088 || $1 == 10089 || $1 == 10090 || $1 == 10091 || $1 == 10091 || $1 == 10092 || $1 == 10093) print;}' > creating_citation_counts_tsv/data/mouse/gene2pubmed
# Getting Rat gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 10114 || $1 == 10115 || $1 == 10116 || $1 == 10117 || $1 == 10118 || $1 == 10119 || $1 == 10121 || $1 == 10122 || $1 == 10127) print;}' > creating_citation_counts_tsv/data/rat/gene2pubmed
# ADD Extra Species gene2pubmed

# Remove main file to save space
rm creating_citation_counts_tsv/data/gene2pubmed

# Next we are going to get refences files for different species :https://hgdownload.soe.ucsc.edu/downloads.html
# To get a list of genomic positions from UCSC to plot the locations of the genes.
# Getting Human refences file
wget -N -P creating_citation_counts_tsv/data/human https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gunzip creating_citation_counts_tsv/data/human/hg38.refGene.gtf.gz
# Getting Mouse refences file
wget -N -P creating_citation_counts_tsv/data/mouse https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz
gunzip creating_citation_counts_tsv/data/mouse/refGene.gtf.gz 
# Getting Rat refences file
wget -N -P creating_citation_counts_tsv/data/rat https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/genes/rn6.refGene.gtf.gz
gunzip creating_citation_counts_tsv/data/rat/rn6.refGene.gtf.gz
# ADD Extra Species refences file

# We need a mapping of gene IDs, which are just numbers, to more meaningful names and descriptions.
wget -N -P creating_citation_counts_tsv/data ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
gunzip creating_citation_counts_tsv/data/gene_info.gz

# Gene info is a big file, so we are parsing out what we need
# Getting Human gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 9606 || $1 == 9605) print;}' > creating_citation_counts_tsv/data/human/gene_info
# Getting Mouse gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 10088 || $1 == 10089 || $1 == 10090 || $1 == 10091 || $1 == 10091 || $1 == 10092) print;}' > creating_citation_counts_tsv/data/mouse/gene_info
# Getting Rat gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 10114 || $1 == 10115 || $1 == 10116 || $1 == 10117 || $1 == 10118 || $1 == 10119 || $1 == 10121 || $1 == 10122 || $1 == 10127) print;}' > creating_citation_counts_tsv/data/rat/gene_info
# ADD Extra Species gene_info

# Remove main file to save space
rm creating_citation_counts_tsv/data/gene_info

# Lastly create the tsv with the total citations per gene along with the gene's information
python3 creating_citation_counts_tsv/scripts/summarize_gene_citations_all_species.py

