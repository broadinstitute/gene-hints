#!/bin/bash

# Creating Citiation TSVs Per selected species. Adapted from source code: https://github.com/pkerpedjiev/gene-citation-counts

# First install 
pip3 install -r creating_citation_counts_tsv/requirements.txt 

# Create a data, tsv, and tmp-ssv folder
mkdir creating_citation_counts_tsv/data
mkdir tsv
mkdir creating_citation_counts_tsv/data/tmp-ssv

#To get an idea of how citation counts have changed over a time frame, we need a list of which citations were published when. 
# This script queries Entrez gene for citation lists for each day for the time frame specified.
days_in_timeframe=180
days_in_timeframe_doubled=$(($days_in_timeframe * 2))
recent_timeframe_end_date=$(date +%Y/%m/%d)
recent_timeframe_start_date=$(date -v -${days_in_timeframe}d +%Y/%m/%d)
past_timeframe_end_date=$(date -v -${days_in_timeframe}d +%Y/%m/%d)
past_timeframe_start_date=$(date -v -${days_in_timeframe_doubled}d +%Y/%m/%d)
output_dir_recent=creating_citation_counts_tsv/data/tmp-ssv/recent-timeframe
output_dir_past=creating_citation_counts_tsv/data/tmp-ssv/past-timeframe

python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --startdate ${recent_timeframe_start_date} --enddate ${recent_timeframe_end_date} --output-dir ${output_dir_recent}
python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --startdate ${past_timeframe_start_date} --enddate ${past_timeframe_end_date} --output-dir ${output_dir_past}

# Consolidate the publications from each day into one complete list
recent_timeframe__year_pmid__ssv_path="creating_citation_counts_tsv/data/recent_pmid_year.ssv"
past_timeframe__year_pmid__ssv_path="creating_citation_counts_tsv/data/past_pmid_year.ssv"
# Cleanup any existing files
rm -f ${recent_timeframe__year_pmid__ssv_path}
rm -f ${past_timeframe__year_pmid__ssv_path}
# For easier processing, we'll collapse the per-day list of files into one file
for file in $(find ${output_dir_recent} -name "*.ssv"); \
do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> ${recent_timeframe__year_pmid__ssv_path}; \
done
for file in $(find ${output_dir_past} -name "*.ssv"); \
do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> ${past_timeframe__year_pmid__ssv_path}; \
done

# Remove tmp-ssv to save space
rm -rf creating_citation_counts_tsv/data/tmp-ssv

# Creating taxonomy folder 
mkdir creating_citation_counts_tsv/taxonomy

# Creating Species folders (human, mouse, rat)
mkdir creating_citation_counts_tsv/data/human
mkdir creating_citation_counts_tsv/data/mouse
mkdir creating_citation_counts_tsv/data/rat
mkdir creating_citation_counts_tsv/data/dog
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
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 9606) print;}' > creating_citation_counts_tsv/data/human/gene2pubmed
# Getting Mouse gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 10090) print;}' > creating_citation_counts_tsv/data/mouse/gene2pubmed
# Getting Rat gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 10116) print;}' > creating_citation_counts_tsv/data/rat/gene2pubmed
# Getting dog gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 9615) print;}' > creating_citation_counts_tsv/data/dog/gene2pubmed
# ADD Extra Species gene2pubmed

# Remove main file to save space
rm creating_citation_counts_tsv/data/gene2pubmed

# Next we are going to get refences files for different species :https://hgdownload.soe.ucsc.edu/downloads.html
# To get a list of genomic positions from UCSC to plot the locations of the genes.
# Getting Human refences file
wget -N -P creating_citation_counts_tsv/data/human https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/human/hg38.refGene.gtf.gz
# Getting Mouse refences file
wget -N -P creating_citation_counts_tsv/data/mouse https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/mouse/refGene.gtf.gz 
# Getting Rat refences file
wget -N -P creating_citation_counts_tsv/data/rat https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/genes/rn6.refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/rat/rn6.refGene.gtf.gz
# Getting dog refences file
wget -N -P creating_citation_counts_tsv/data/dog https://hgdownload.soe.ucsc.edu/goldenPath/canFam5/bigZips/genes/refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/dog/refGene.gtf.gz
# ADD Extra Species refences file

# We need a mapping of gene IDs, which are just numbers, to more meaningful names and descriptions.
wget -N -P creating_citation_counts_tsv/data ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
gunzip creating_citation_counts_tsv/data/gene_info.gz

# Gene info is a big file, so we are parsing out what we need
# Getting Human gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 9606) print;}' > creating_citation_counts_tsv/data/human/gene_info
# Getting Mouse gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 10090) print;}' > creating_citation_counts_tsv/data/mouse/gene_info
# Getting Rat gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 10116) print;}' > creating_citation_counts_tsv/data/rat/gene_info
# Getting dog gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 9615) print;}' > creating_citation_counts_tsv/data/dog/gene_info
# ADD Extra Species gene_info

# Remove main file to save space
rm creating_citation_counts_tsv/data/gene_info

# Lastly create the tsv with the total citations per gene along with the gene's information
python3 creating_citation_counts_tsv/scripts/summarize_gene_citations_all_species.py ${recent_timeframe__year_pmid__ssv_path} ${past_timeframe__year_pmid__ssv_path} $days_in_timeframe
