

import os
import csv
from datetime import date, timedelta, datetime

today = date.today().strftime("%Y/%m/%d")
month_ago = (datetime.today() - timedelta(days=155)).strftime("%Y/%m/%d")

f = open("data/ucsc-data/refGene.txt", "r")
list_append = {}
for lines in f.readlines():
    # if 'transcript' in lines.split("\t")[2]:
        # gene_name = (lines.split("\t")[8].split(";")[0].replace("gene_id ", "").replace('"', '')).upper()
        # chr_id = lines.split("\t")[0].replace("chr","")
        # start = lines.split("\t")[3]
        # lenght = abs(int(lines.split("\t")[3])-int(lines.split("\t")[4]))
        # list_append[gene_name] = [gene_name, chr_id, start, lenght]
    list_append[(lines.split("\t")[12]).upper()] = [(lines.split("\t")[12]).upper(), lines.split("\t")[2].replace("chr",""), lines.split("\t")[4],abs(int(lines.split("\t")[4])-int(lines.split("\t")[5])) ]


directory = 'all_species_tsvs'
name = []
top_ten_citation = {}
# os.listdir(directory)
for filename in ['homo-sapiens-citation-information.tsv']:
    if filename.endswith(".tsv"):
        new_rows = []
        tsv_file = open(os.path.join(directory, filename), encoding = "ISO-8859-1")
        read_tsv = csv.reader(tsv_file, delimiter="\t")
        count = 0

        for row in read_tsv:
            # print(row[0].upper())
            result = [x for x in list_append.keys() if x.startswith(row[0].upper())] 
            # print(result)
            if "#name" in row[0]:
                new_rows.append(row)
            if result:
                row[1] = list_append[result[0]][1]
                row[2] = list_append[result[0]][2]
                row[3] = list_append[result[0]][3]
                new_rows.append(row)
        with open(filename, 'wt') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            for new_row in new_rows:
                tsv_writer.writerow(new_row)
            
        top_ten_citation[filename.replace("_citation_information.tsv","")] = new_rows
exit(1)
with open('top_ten_citation_per_species.html', 'w') as f:
    html_text = "<h1> The Top Ten Citation Per Species</h1><br>"
    for tile_name in top_ten_citation.keys():
        title = tile_name.replace("_", " ").upper()
        html_text = html_text + f"<h2>{title}</h2> <table style='width:100%'>"
        html_text = html_text + "<tr> <th>name</th> <th>chromosome</th> <th>start</th> <th>length</th> <th>color</th> <th>full_name</th> <th>Citations</th> <th>signifance</th></tr>"
        for line in top_ten_citation[tile_name]:
            html_text = html_text + "<tr>"
            for column in line:
                html_text = html_text + f"<td>{column}</td>"
            html_text = html_text + "</tr>"
        html_text = html_text + "</table>"
    f.write(html_text)

