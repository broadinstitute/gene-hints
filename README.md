# gene-hints
Discoverability for gene search

[Gene Hints](https://broad.io/gene-hints) is a reusable UI component backed by genomic data.  It makes gene search more discoverable by finding and showing interesting genes across the genome.  Users can then click genes to search them.

Gene Hints complements existing gene search in any web application, and requires no backend to operate.

# User interface
## Show interesting genes by default
A curated set of human gene hints shows genes of interest to search, letting users discover genomic landmarks at a glance.
![Gene Hints](https://raw.githubusercontent.com/broadinstitute/gene-hints/main/img/01-gene-hints.png)

## Tooltips give more context
Hover over the ACE2 gene to see its significance.  COVID-19 virus enters cells through the ACE2 protein.
![Gene Hints, ACE2 tooltip](https://raw.githubusercontent.com/broadinstitute/gene-hints/main/img/02-gene-hints-ace2-tooltip.png)

## Hints for downstream exploration
Click ACE2 to search it, then see similar or interacting genes.
![Gene Hints, ACE2 related genes](https://raw.githubusercontent.com/broadinstitute/gene-hints/main/img/03-gene-hints-ace2-related-genes.png)

## Adapt to new genomic context
A tooltip for a related gene, AGT, shows the pathway in which it interacts with the searched gene.
![Gene Hints, ACE2 relate gene AGT tooltip](https://raw.githubusercontent.com/broadinstitute/gene-hints/main/img/04-gene-hints-ace2-related-genes-agt-tooltip.png)

## Support for multiple species
Gene Hints also works for mouse and rat.  A tooltip shows significance and citation counts for Lepr, a popular mouse gene.
![Gene Hints, mouse LEPR gene tooltip](https://raw.githubusercontent.com/broadinstitute/gene-hints/main/img/05-gene-hints-mus-musculus-lepr-tooltip.png)
