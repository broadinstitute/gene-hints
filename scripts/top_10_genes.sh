#!/bin/bash
#TODO: ADD script README and install.py
# Get NCBI's genes to publications table, source code: https://github.com/pkerpedjiev/gene-citation-counts
mkdir -p data/tmp-ssv
mkdir -p data/genbank-data
mkdir -p data/ucsc-data
mkdir -p data/taxonomy
mkdir -p all_species_tsv

#TODO install csv
#To get an idea of how citation counts have changed over the years, 
#we need a list of which citations were published when. 
# This script queries Entrez gene for citation lists for each day.
python3 scripts/pmids_by_date.py

#Consolidate the publications from each day into one complete list
#For easier processing, we'll collapse the per-day list of files into one file
for file in $(find . -name "*.ssv"); \
do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> data/genbank-data/recent_pmid_year.ssv; \
done

rm data/tmp-ssv/*
rmdir data/tmp-ssv

# The NCBI maintains a list of citations that mention each gene in its database. We'll download it.
# TODO: Download wget
wget -N -P data/genbank-data/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz
gunzip data/genbank-data/gene2pubmed.gz
mkdir -p data/ucsc-data
 
# We'll get a list of chromosome sizes from the UCSC genome browser for plotting purposes:
wget -N -P data/ucsc-data/ http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz

# We'll also grab a list of genomic positions from UCSC to plot the locations of the genes.
wget -N -P data/ucsc-data/ http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz

# We need a mapping of gene IDs, which are just numbers, to more meaningful names and descriptions.
wget -N -P data/genbank-data/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
gunzip data/genbank-data/gene_info.gz

wget -N -P data/taxonomy/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z 
zcat data/taxonomy/new_taxdump.tar.Z > data/taxonomy/new_taxdump.tar; /
tar -C data/taxonomy/ -xvf data/taxonomy/new_taxdump.tar
mv data/taxonomy/names.dmp data/genbank-data/taxonomy
rm data/taxonomy/*


# Get gene to refseq information so that we can use it to get transcript level information. 
# NCBI's gene database annotates only genes. Actual transcript information is stored in RefSeq.
wget -N -P data/genbank-data/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
gunzip data/genbank-data/gene2refseq.gz

# All of the aggregations which count 
# the total citations per gene as well as the citations 
# per per gene per year are calculated using the summarize_gene_citations.py script. 
# It could be done using command line scripts but I started with python so I'll stick with it.
# aggregate by citation counts and by year
# TODO install pyspark and shortuuid https://www.apache.org/dyn/closer.lua/spark/spark-3.0.1/spark-3.0.1-bin-hadoop2.7.tgz     https://stackoverflow.com/questions/30518362/how-do-i-set-the-drivers-python-version-in-spark
# instal install psutil 
python3 scripts/summarize_gene_citations_all_species.py

# rm data/genbank-data/*
# rm data/ucsc-data/*
# rmdir data/ucsc-data
# rmdir data/genbank-data