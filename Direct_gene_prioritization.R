
#Load libraries

if (!require("data.table")) install.packages("data.table", dependencies = T)
if (!require("tidyverse")) install.packages("tidyverse", dependencies = T)
if (!require("ontologyIndex")) install.packages("ontologyIndex", dependencies = T)
if (!require("ontologySimilarity")) install.packages("ontologySimilarity", dependencies = T)


args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments (TXT filename and HPO terms) are provided
if (length(args) < 3) {
  cat("Usage: Rscript direct_gene_prioritization.R <TXT_filename> <HPO_terms> <Output_filename>\n")
  quit(status = 1)
}

df_hpo <- read_tsv("hpo2023.txt", show_col_types = FALSE)

# Read ontology file 
hpo <-  get_ontology("hpo2023.obo")

#Create HPO gene list
hpo_genes <- split(df_hpo$`HPO-Term-ID`, df_hpo$Gene_symbol)

txt_filename <- args[1]
output_filename <- args[length(args)]
query_set <- args[-c(1, length(args))]

# Read the CSV file into a data frame
tab <- read.table(txt_filename, header=FALSE, col.names = "GeneName", as.is = TRUE)

tab$GeneName <- toupper(tab$GeneName)
colnames(tab) <- "Gene"
# Take HP terms from user and Split query set
query_set <- unlist(strsplit(str_replace_all(query_set, " ", ""),split=','))

sim_grid <- get_sim_grid(ontology=hpo, term_sets=hpo_genes, term_sets2=list(query_set), combine = "average")

print("Calculating scores for each gene, please wait...")

df.hpo <- as.data.frame(sim_grid)
df.hpo <- tibble::rownames_to_column(df.hpo, "Gene")
colnames(df.hpo)[2] <- "OntologySimilarity_Score"
df.hpo_all <- df.hpo %>% arrange(desc(df.hpo$OntologySimilarity_Score))
df.hpo_all$OntologySimilarity_Score <- round(df.hpo_all$OntologySimilarity_Score,5)

dt <- left_join(tab, df.hpo_all)

write.csv(dt, file = output_filename, row.names = F)






  
  
