####################################################################################################################################################
#This script was written by Christina Canavati - the Tabach laboratory 
#This script receieves patient HPO terms and a list of genes as input and outputs:
#A ranked gene list according to phylogenetic profiling similarity and STRING interaction with HPO-related genes (retrieved from a prebuilt dataset)
####################################################################################################################################################

#Load libraries

if (!require("data.table")) install.packages("data.table", dependencies = T)
if (!require("tidyverse")) install.packages("data.table", dependencies = T)
if (!require("ontologyIndex")) install.packages("data.table", dependencies = T)
if (!require("ontologySimilarity")) install.packages("data.table", dependencies = T)
if (!require("survcomp")) install.packages("data.table", dependencies = T)
if (!require("readr")) install.packages("data.table", dependencies = T)

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments (TXT filename and HPO terms) are provided
if (length(args) < 3) {
  cat("Usage: Rscript gene_prioritization.R <TXT_filename> <HPO_terms> <Output_filename>\n")
  quit(status = 1)
}


df_hpo <- read_tsv("hpo2022.txt", show_col_types = FALSE)

# Read ontology file 
hpo <-  get_ontology("hpo2022.obo")

#Create HPO gene list
hpo_genes <- split(df_hpo$`HPO-Term-ID`, df_hpo$Gene_symbol)

#Read All_Clades_File and clades index files
All_Clades_File <- fread("gunzip -c UNIPROT_ALL_CLADES_FILE.csv.gz")
ind <- fread("Clade_Indeces_File.csv")
string <- fread("gunzip -c STRING_Database_HPO.csv.gz")

#Read user inquery
# Get the CSV filename and HPO terms from command-line arguments
txt_filename <- args[1]
output_filename <- args[length(args)]
query_set <- args[-c(1, length(args))]

# Read the CSV file into a data frame
tab <- read.table(txt_filename, header=FALSE, col.names = "GeneName", as.is = TRUE)

  genes <- tab$GeneName
  genes <- toupper(genes)
  genes <- unlist(strsplit(genes, ",", " "))

  # Take HP terms from user and Split query set
  query_set <- unlist(strsplit(str_replace_all(query_set, " ", ""),split=','))
  
  sim_grid <- get_sim_grid(ontology=hpo, term_sets=hpo_genes, term_sets2=list(query_set), combine = "average")
  
  # Convert matrix into a dataframe                 
  df.hpo <- as.data.frame(sim_grid)
  df.hpo <- tibble::rownames_to_column(df.hpo, "Gene")
  colnames(df.hpo)[2] <- "Score"
  df.hpo_all <- df.hpo %>% arrange(desc(df.hpo$Score))
  df.hpo_all$Score <- round(df.hpo_all$Score,5)
  colnames(df.hpo_all)[1] <- "Gene2"
  
  dt <- data.frame(Gene = 0, p_value = 0)#, sample_ID = 0)
  dt <- dt[NULL,]

print("Calculating scores for each gene, please wait...")

suppressMessages({
suppressWarnings({  

  for(gene in genes){
   
    Gene_Rank_File <- All_Clades_File[Gene1 == gene]
    Gene_Rank_File <- left_join(Gene_Rank_File, ind, by ="Clade_Index")
    Gene_Rank_File <- left_join(Gene_Rank_File, df.hpo_all)
    Gene_Rank_File <- Gene_Rank_File[Rank <= 50] %>% mutate(Clade = "Coevolution")
    Gene_Rank_File <- Gene_Rank_File[,-c(1,4)]
    
    string_gene <- string[P1 == gene]
    colnames(string_gene)[3] <- "Gene2"
    string_gene <- inner_join(string_gene, df.hpo_all)
    string_gene <- data.table(string_gene)
    string_gene <- string_gene[,-1]
    colnames(string_gene) <- c("Gene1", "Gene2", "Score")
    
    if(nrow(string_gene) > 0){
      string_gene <- string_gene %>% mutate(Clade = "STRING")
      Gene_Rank_File <- rbind(Gene_Rank_File, string_gene)}
    
    savethis <- data.frame(Gene = 0, Clade = 0, p_value =0)
    savethis <- savethis[NULL,]

    for(i in unique(Gene_Rank_File$Clade)){
      
      clade <- Gene_Rank_File[Clade == i]
      if(nrow(clade) == 0){
        next
      }

      test <- ks.test(clade$Score, df.hpo_all$Score, alternative = "less")
      p_value = test$p.value
      tmp <- data.frame(Gene = gene, Clade = i, p_value = p_value)
      
      savethis <- rbind(savethis, tmp)
    }
    
    combined_p_value = combine.test(savethis$p_value, method = "fisher")
    if(length(combined_p_value) == 0){combined_p_value = 1}
    dt_tmp <- data.frame(Gene = gene, p_value = combined_p_value)#, sample_ID = tools::file_path_sans_ext(n))
    dt <- rbind(dt, dt_tmp)
  }
  })
})

print(paste0("Creating output file ", output_filename))

  dt <- dt %>% mutate(Rank = frank(p_value, ties.method = "dense"))
  colnames(dt) <- c("Gene", "EvORanker p-value", "Rank")
  write.csv(dt, file = output_filename, row.names = F)
  
