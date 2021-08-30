#!/bin/bash

# Creates citiation TSVs for selected organisms
# Adapted from https://github.com/pkerpedjiev/gene-citation-counts

# First install
pip3 install -r creating_citation_counts_tsv/requirements.txt

# Create a data, tsv, and tmp folder
mkdir creating_citation_counts_tsv/data
mkdir tsv
mkdir creating_citation_counts_tsv/data/tmp

# To get an idea of how citation counts have changed over a time frame, we need
# a list of which citations were published when.  This script queries the NCBI
# Entrez Gene database for citation lists for each day for the time frame
# specified.
days_in_timeframe=60
days_in_timeframe_doubled=$(($days_in_timeframe * 2))
timeframe_end_date=$(date +%Y/%m/%d)
timeframe_start_date=$(date -v -${days_in_timeframe}d +%Y/%m/%d)
prev_timeframe_end_date=$(date -v -${days_in_timeframe}d +%Y/%m/%d)
prev_timeframe_start_date=$(date -v -${days_in_timeframe_doubled}d +%Y/%m/%d)
output_dir=creating_citation_counts_tsv/data/tmp/recent-timeframe
prev_output_dir=creating_citation_counts_tsv/data/tmp/past-timeframe

python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --start-date ${timeframe_start_date} --end-date ${timeframe_end_date} --output-dir ${output_dir}
python3 creating_citation_counts_tsv/scripts/pmids_by_date.py --start-date ${prev_timeframe_start_date} --end-date ${prev_timeframe_end_date} --output-dir ${prev_output_dir}

# Consolidate the publications from each day into one complete list
pmid_times_path="creating_citation_counts_tsv/data/pmid_times.tsv"
prev_pmid_times_path="creating_citation_counts_tsv/data/prev_pmid_times.tsv"
# Cleanup any existing files
rm -f ${pmid_times_path}
rm -f ${prev_pmid_times_path}
# For easier processing, we'll collapse the per-day list of files into one file
for file in $(find ${output_dir} -name "*.tsv"); \
do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> ${pmid_times_path}; \
done
for file in $(find ${prev_output_dir} -name "*.tsv"); \
do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> ${prev_pmid_times_path}; \
done

# Remove tmp to save space
rm -rf creating_citation_counts_tsv/data/tmp

# Create taxonomy folder
mkdir creating_citation_counts_tsv/taxonomy

# Create organism folders
mkdir creating_citation_counts_tsv/data/human
mkdir creating_citation_counts_tsv/data/mouse
mkdir creating_citation_counts_tsv/data/rat
mkdir creating_citation_counts_tsv/data/dog
mkdir creating_citation_counts_tsv/data/cat

# Get scientific organism names
wget -N -P creating_citation_counts_tsv/taxonomy/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z
zcat creating_citation_counts_tsv/taxonomy/new_taxdump.tar.Z > creating_citation_counts_tsv/taxonomy/new_taxdump.tar; /
tar -C creating_citation_counts_tsv/taxonomy/ -xvf creating_citation_counts_tsv/taxonomy/new_taxdump.tar
mv creating_citation_counts_tsv/taxonomy/names.dmp creating_citation_counts_tsv/taxonomy_name
rm creating_citation_counts_tsv/taxonomy/*

# NCBI maintains a list of citations that mention each gene in its database.  We'll download it.
wget -N -P creating_citation_counts_tsv/data ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz
gunzip creating_citation_counts_tsv/data/gene2pubmed.gz

# gene2pubmed is a big file, so we are parsing out what we need
# Get human gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 9606) print;}' > creating_citation_counts_tsv/data/human/gene2pubmed
# Get mouse gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 10090) print;}' > creating_citation_counts_tsv/data/mouse/gene2pubmed
# Get rat gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 10116) print;}' > creating_citation_counts_tsv/data/rat/gene2pubmed
# Get dog gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 9615) print;}' > creating_citation_counts_tsv/data/dog/gene2pubmed
# Get cat gene2pubmed
cat creating_citation_counts_tsv/data/gene2pubmed | awk '{if ($1 == 9685) print;}' > creating_citation_counts_tsv/data/cat/gene2pubmed

# Remove main file to save space
rm creating_citation_counts_tsv/data/gene2pubmed

# Next we'll get reference files for different species: https://hgdownload.soe.ucsc.edu/downloads.html
# to get a list of genomic positions from UCSC to plot the locations of the genes.
# Get human references file
wget -N -P creating_citation_counts_tsv/data/human https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/human/hg38.refGene.gtf.gz
# Get human references file
wget -N -P creating_citation_counts_tsv/data/mouse https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/mouse/refGene.gtf.gz
# Get human references file
wget -N -P creating_citation_counts_tsv/data/rat https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/genes/rn6.refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/rat/rn6.refGene.gtf.gz
# Get dog references file
wget -N -P creating_citation_counts_tsv/data/dog https://hgdownload.soe.ucsc.edu/goldenPath/canFam5/bigZips/genes/refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/dog/refGene.gtf.gz
# Get cat references file
wget -N -P creating_citation_counts_tsv/data/cat https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/genes/felCat9.refGene.gtf.gz
gunzip -f creating_citation_counts_tsv/data/cat/felCat9.refGene.gtf.gz

# We need a map of gene IDs, which are just numbers, to more meaningful names and descriptions.
wget -N -P creating_citation_counts_tsv/data ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
gunzip creating_citation_counts_tsv/data/gene_info.gz

# Gene info is a big file, so we are parsing out what we need
# Get human gene_info
cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 9606) print;}' > creating_citation_counts_tsv/data/human/gene_info
# # Get mouse gene_info
# cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 10090) print;}' > creating_citation_counts_tsv/data/mouse/gene_info
# # Get rat gene_info
# cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 10116) print;}' > creating_citation_counts_tsv/data/rat/gene_info
# # Get dog gene_info
# cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 9615) print;}' > creating_citation_counts_tsv/data/dog/gene_info
# # Get cat gene_info
# cat creating_citation_counts_tsv/data/gene_info | awk '{if ($1 == 9685) print;}' > creating_citation_counts_tsv/data/cat/gene_info

# Remove main file to save space
rm creating_citation_counts_tsv/data/gene_info

# Lastly create the TSV with the total citations per gene along with the gene's information

echo "python3 creating_citation_counts_tsv/scripts/summarize_gene_citations_all_species.py ${pmid_times_path} ${prev_pmid_times_path} $days_in_timeframe"
python3 creating_citation_counts_tsv/scripts/summarize_gene_citations_all_species.py ${pmid_times_path} ${prev_pmid_times_path} $days_in_timeframe
