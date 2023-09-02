# EvORanker
EvORanker is a software that implements phylogenetic profiling and STRING-based prioritization of candidate  genes for human diseases.

## Introduction
EvORanker is the first platform that leverages a clade-wise phylogenetic profiling (PP) approach to associate gene variants in patients with their phenotypes. It is particularly useful for identifying disease associations for genes that may not have been previously characterized or published.

## Releases
EvORanker is also available as a web-sever at [here](https://ccanavati.shinyapps.io/EvORanker/).

## Dependencies

- R (version 3.6 or higher)

## Installation

Please clone the repository into your computer:
```
git clone https://github.com/WGLab/phenolyzer](https://github.com/ccanavati/EvoRanker
```
Then enter the EvoRanker directory:
```
cd EvoRanker
```
## Synopsis

- Prioritize genes based on PP and STRING (The smaller the p-value, the stronger the evidence against the null hypothesis, indicating a more significant association for the gene): 
```
Rscript gene_prioritization.R example_genes.txt HP:0001188 HP:0001166 HP:0001199 output_file.csv
```

## Contact
- Christina Canavati (canavatichristina@gmail.com)
