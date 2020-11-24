import os
import findspark
import pyspark
import shortuuid
import shutil
import time
import os
import csv

# 
findspark.init()

# 
sc = pyspark.SparkContext()

#
shortuuid.uuid()

#
assembly = ''

#
op = os.path

#
data_dir = op.expanduser("data")
output_dir = op.join(data_dir, assembly)   # where all of the intermediate output will be stored
base_ucsc_dir = op.join(data_dir, 'ucsc-data/{}'.format(assembly))  # where all of the files downloaded from UCSC will be stored



# create a directory to store intermediate output files
def get_outfile(table_name):
    outfile = op.join(output_dir, 'genbank-output/{}'.format(table_name))
    if op.exists(outfile):
        shutil.rmtree(outfile)
    return outfile

# Loading the refgene data
base_dir=op.join(op.expanduser("data/genbank-data/"), assembly)

gene2pubmed = (sc.textFile(op.join(base_dir, "gene2pubmed"))
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .map(lambda x: ((int(x[0]), int(x[1])), int(x[2])))
                     )
print(gene2pubmed.take(1))

taxonomy = (sc.textFile(op.join(base_dir, "taxonomy"))
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .map(lambda x: ((str(x[0]), str(x[2]), str(x[6]))))
                     )
print(taxonomy.take(1))

ucsc_refGene = (sc.textFile(op.join(base_ucsc_dir, "refGene.txt"))
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .map(lambda x: ((str(x[2]).replace("chr",""), abs(int(x[4]) - int(x[5])), str(x[12]))))
                     )
print(ucsc_refGene.take(1))

ucsc_refGene_list = ucsc_refGene.take(ucsc_refGene.count())
ucsc_refGene_dict = {}
for ucsc_refGene_pair in ucsc_refGene_list:
    if ucsc_refGene_pair[2].upper() not in ucsc_refGene_dict.keys():
        ucsc_refGene_dict[ucsc_refGene_pair[2].upper()] = [ucsc_refGene_pair[0], ucsc_refGene_pair[1]]

taxonomy_list = taxonomy.take(taxonomy.count())
taxonomy_dict = {}
for taxonomy_pair in taxonomy_list:
    if taxonomy_pair[0] not in taxonomy_dict.keys():
        taxonomy_dict[taxonomy_pair[0]] = taxonomy_pair[1]
    elif "scientific name" in taxonomy_pair[2]:
        taxonomy_dict[taxonomy_pair[0]] = taxonomy_pair[1]

t1 = time.time()
pubmeds_set = set([x[1] for x in gene2pubmed.collect()])
t2 = time.time()
print("time taken", t2 - t1)

from pyspark import SparkContext, SparkConf

base_dir = op.join(data_dir, 'genbank-data/')
taxid_gene_info = (sc.textFile(op.join(base_dir, 'gene_info'))
                   .filter(lambda x: x[0] != '#')
                   .map(lambda x: x.split('\t'))
                   .map(lambda x: ((int(x[0]), int(x[1])),(x[2], x[8], x[9])))
                   )
print(taxid_gene_info.take(1))

taxid_gene_refseq_id = (sc.textFile(op.join(base_dir, "gene2refseq"))
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .map(lambda x: ((int(x[0]), int(x[1]), str(x[9]), str(x[10])), (x[3].split('.')[0])))
                        )
print(taxid_gene_refseq_id.take(1))

taxid_gene_refseq_id_list = taxid_gene_refseq_id.take(taxid_gene_refseq_id.count())
taxid_gene_refseq_id_dict = {}
for taxid_gene_refseq_id_pair in taxid_gene_refseq_id_list:
    if taxid_gene_refseq_id_pair[0][1] not in taxid_gene_refseq_id_dict.keys():
        taxid_gene_refseq_id_dict[str(taxid_gene_refseq_id_pair[0][1])] = [taxid_gene_refseq_id_pair[0][2], taxid_gene_refseq_id_pair[0][3]]
import time
t1 = time.time()
taxid_gene_info_refseq = taxid_gene_info.join(taxid_gene_refseq_id)

t2 = time.time()
print("time taken", t2 - t1)

len(pubmeds_set)
year_pmid = (sc.textFile(op.join(base_dir, 'recent_pmid_year.ssv'))
                  .map(lambda x: x.split())
                  .map(lambda x: (int(x[0]), int(x[1])))
                  .filter(lambda x: x[1] in pubmeds_set))
year_pmid_collected = year_pmid.collect()

pmid_year = dict([(x[1], x[0]) for x in year_pmid_collected])

print([k for k in list(pmid_year.values())[:10]])

taxid_gene_info_pubmed = (taxid_gene_info.join(gene2pubmed)
                    .filter(lambda x: x[1][1] in pmid_year)
                                 .map(lambda x: ((x[0][0], x[0][1], x[1][1]), x[1])))
print(taxid_gene_info_pubmed.take(1))        

# get total citation counts over all time for each taxid, gene_id combo
gene_counts = sorted(taxid_gene_info_pubmed.map(lambda x: ((x[0][0], x[0][1]), (x[1][0], 1)))
.reduceByKey(lambda x1, x2: (x1[0], x1[1] + x2[1]))
.collect(), key=lambda x: -x[1][1])

gene_counts[:3]
print(gene_counts[:3])

taxid_gene_info_year = (taxid_gene_info_pubmed.map(lambda x: ((x[0][0], x[0][1], pmid_year[x[1][1]]), x[1]))
                        .map(lambda x: ((x[0][0], x[0][1], pmid_year[x[1][1]]), (x[1][0], 1)))
                                 .reduceByKey(lambda x1, x2: (x1[0], x1[1] + x2[1])))
print(taxid_gene_info_year.take(1))
print(taxid_gene_info_year.count())
print(taxid_gene_info_year.filter(lambda x: x[1][0][2] == 'rRNA').take(1))

genes_per_year = dict(taxid_gene_info_pubmed.map(lambda x: ((pmid_year[x[1][1]], x[0][1]), 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .map(lambda x: (x[0][0], 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .collect())

citations_per_year = dict((taxid_gene_info_pubmed.map(lambda x: (pmid_year[x[1][1]], 1))
                      .reduceByKey(lambda x1,x2: x1+x2)
                      .collect()))

citation_genes = sorted(taxid_gene_info_pubmed.map(lambda x: (x[1][1], 1))
 .reduceByKey(lambda x1, x2: x1 + x2)
 .collect(), key=lambda x: -x[1])

print("citation_genes:", citation_genes[:10])

gene_types_citations = (taxid_gene_info_pubmed.map(lambda x: (x[1][0][2], 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .collect()
)
print("gene_types_citations:", gene_types_citations)

gene_types_counts = (taxid_gene_info_pubmed.map(lambda x: ((x[1][0][2], x[0][1]),1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .map(lambda x: (x[0][0], 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .collect())

for gtc in gene_types_counts:
    print(gtc[0], gtc[1])
string_values = (taxid_gene_info_year.map(lambda x: "\t".join(map(str, 
                 [x[0][0], x[0][1], x[0][2], x[1][0][0], x[1][0][1], x[1][0][2], x[1][1]
                 ])))
                 .collect())
gene_id_dict = {}

for values in string_values:
    values_list = values.split("\t")
    citation = values_list[-1]
    gene_id = values_list[0]
    gene_name = values_list[3]
    gene_full_name = values_list[4]
    gene_chromosome = "N/a"
    gene_length = "N/A"
    gene_start = "N/A"
    if gene_name.upper() in ucsc_refGene_dict.keys():
        print("in")
        gene_chromosome = ucsc_refGene_dict[gene_name][0]
        gene_length = ucsc_refGene_dict[gene_name][1]
    if gene_id in taxid_gene_refseq_id_dict.keys(): 
        print("in in")
        gene_start = taxid_gene_refseq_id_dict[gene_id][0]
        gene_length = taxid_gene_refseq_id_dict[gene_id][1]
        print(gene_start)
    
    new_values = f"{gene_name}\t{gene_chromosome}\t{gene_start}\t{gene_length}\t#73af42\t{gene_full_name}\t{citation}"
    if gene_id not in gene_id_dict.keys():
        gene_id_dict[gene_id] = [new_values]
    else: 
        gene_id_dict[gene_id].append(new_values)

top_ten_citation = {}
for gene_id_tax_id in gene_id_dict.keys():
    rows = []
    if gene_id_tax_id in taxonomy_dict.keys():
        gene_species_name = taxonomy_dict[gene_id_tax_id].replace(" ", "_").replace("(", "_").replace(")", "_").replace("/", "_").replace("=", "_")
        tsv_name= f'all_species_tsv/{gene_species_name}_citation_information.tsv'
        print(tsv_name)    
        for info in gene_id_dict[gene_id_tax_id]:
            info_list = info.split("\t")
            rows.append(info_list)
        rows.sort(key=lambda x: x[-1])
        top_ten_rows = rows[-10:][::-1]
        with open(tsv_name, 'wt') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow(["#name", "chromosome", "start", "length", "color", "full_name"])
            for ttr in top_ten_rows[:-1]:
                tsv_writer.writerow(ttr)
        top_ten_citation[gene_species_name] = top_ten_rows
    print("gene_id_tax_id")
    print(gene_id_tax_id)


with open('top_ten_citation_per_species.html', 'w') as f:
    html_text = "<h1> The Top Ten Citation Per Species</h1><br>"
    for tile_name in top_ten_citation.keys():
        title = tile_name.replace("_", " ").upper()
        html_text = html_text + f"<h2>{title}</h2> <table style='width:100%'>"
        html_text = html_text + "<tr> <th>name</th> <th>chromosome</th> <th>start</th> <th>length</th> <th>color</th> <th>full_name</th> <th>Citations</th></tr>"
        for line in top_ten_citation[tile_name]:
            html_text = html_text + "<tr>"
            for column in line:
                html_text = html_text + f"<td>{column}</td>"
            html_text = html_text + "</tr>"
        html_text = html_text + "</table>"
    f.write(html_text)


