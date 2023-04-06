###############################################################################
#This script was written by Christina Canavati of the Tabach laboratory 
#EvORanker was built using the R shiny package
#This script receieves patient HPO terms and a list of genes as input and ouputs:
#(1) A ranked gene list according to human phenotype similarity using the OntologySimilarity package
#(2) A ranked gene list according to phylogenetic profiling similarity and STRING interaction with HPO-related genes (retrieved from a prebuilt dataset)
#(4) A barplot displaying the scored gene list results
#(5) A coevolution subnetwork showing the HPO-related genes coevolving and STRING interacting with the query gene
#(6) An EnrichR site direction to show the enrichment of pathways and ontologies of the co-evolving and STRING-interacting disease genes
#This script was last updated on September 29, 2022

# Import libraries
library(shiny)
library(shinythemes)
library(data.table)
library(RCurl)
library(shinyBS)
library(tidyverse)
library(ggthemes)
library(ontologyIndex)
library(ontologySimilarity)
library(rjson)
library(shinydashboard)
library(DT)
library(shinyjs)
library(hrbrthemes)
library(plotly)
library(viridis)
library(dashboardthemes)
library(visNetwork)
library(R.utils)
library(httr)
library(BiocManager)
options(repos = BiocManager::repositories())
library(survcomp)
#Read HPO txt file
df_hpo <- read_tsv("hpo2022.txt", show_col_types = FALSE)
colnames(df_hpo) <- c("Gene_ID", "Gene_symbol", "HPO-Term-ID", "HPO-Term-Name", "Frequency-raw", "Frequency-HPO", "Info_Source", "GD_Souce_ID", "Disease ID")

# Read ontology file 
hpo <-  get_ontology("hpo2022.obo")

#Create HPO gene list
hpo_genes <- split(df_hpo$`HPO-Term-ID`, df_hpo$Gene_symbol)

#Read All_Clades_File and index 
All_Clades_File <- fread("gunzip -c UNIPROT_ALL_CLADES_FILE.csv.gz")
ind <- fread("Clade_Indeces_File.csv")

#Columns displayed on DT table8
display_colnames <- c("Gene", "EvORanker p-value")
display_colnames_HPO <- c("Gene", "Combined Semantic Simimlarity Score")
#All NPP genes
All_genes <- unique(All_Clades_File$Gene1)

#Read STRING database
string <- fread("gunzip -c STRING_Database_HPO.csv.gz")

#Read precalculated simulations file
#simulations <- readRDS("Simulations_Full.rds")
IC <- get_term_info_content(hpo, hpo$id)

#Specify empty list names
# plist <- list()
# glist <- list()
# clist <- list()
# hpolist <- list()

#ENRICHR
ENRICHR_ADDLIST_URL <- 'https://maayanlab.cloud/Enrichr/addList'
ENRICHR_DATASET_URL <- "https://maayanlab.cloud/Enrichr/enrich?dataset="

#Color functions
myPal <- colorRampPalette(c("#0868ac","#bae4bc"))
myPal2 <- colorRampPalette(c("#ffffff","#ebdb41", "red"))

js_code <- "
shinyjs.browseURL = function(url) {
  window.open(url,'_blank');
}
"
####################################
# User interface                   #
####################################
ui <- dashboardPage(
  dashboardHeader(title = "EvORanker"),
  dashboardSidebar(collapsed = T, width = 300,
                   useShinyjs(),
                   extendShinyjs(text = js_code, functions = 'browseURL'),
                   ############################
                   htmlOutput("markup"),
                   textOutput("clicked_text"),
                   tags$head(tags$style(HTML('.skin-blue .sidebar-menu > li.active > a,
                                                .skin-blue .sidebar-menu > li:hover > a {
                                                border-left-color: #ff0000;}'))),
                   h4(strong("Coevolution-Based Gene Prioritization"), align = "center"),
                   br(),
                   tags$h4("Insert HPO terms:",style="color:black;text-align:center;"),
                   textInput("txt1", "HPO terms", "", placeholder = "e.g. HP:0001166, HP:0001250"),
                   tags$h4("Insert candidate genes:",style="color:black;text-align:center;"),
                   h6("Use valid HGNC symbols without quotemarks",style="color:blue;text-align:center;font-weight: bold;font-size:14px;"),
                   h6("It is not recommended to exceed 100 genes",style="color:orange;text-align:center;font-size:14px;font-weight: bold;"),
                   textAreaInput("txt2", "Patient Genes", "" , width = "600px", height = "400px", rows = 5, placeholder = "e.g. SCN1A\nCACNA1C\nAPTX"),
                   actionButton("submitbutton", "Submit", class = "btn btn-primary")
  ), #End of dashboardSidebar
  
  dashboardBody(
    
    tags$style(HTML("
.box-header h3.box-title {
   font-weight: bold;
   font-size: 24px;}
.box.box-solid.box-primary>.box-header {
  color:#fff;
  background:#2d67cc}
.box.box-solid.box-primary{
border-bottom-color:#2d67cc;
border-left-color:#2d67cc;
border-right-color:#2d67cc;
border-top-color:#2d67cc;
}
 ")), 
    
    shinyDashboardThemes(
      theme = "grey_light"
    ),
    tabsetPanel(id = "tabs", 
                tabPanel(title = "Home", id="home_tab",
                         fluidRow(style="background:url(https://i.im.ge/2022/08/26/OwFfNh.giphy.gif) no-repeat center;background-size:300px;",
                                  tags$div(class="banner_headline",h2(HTML("<b>EvORanker</b>"), align = "center"))
                         ), #End of fluidRow
                         
                         fluidRow(
                           box(width = 12, 
                               div(align = "center",
                                   h1(HTML("<b>Welcome to EvORanker</b>")), status = "primary", 
                                   h2("EvORanker is an interactive web tool for prioritizing candidate genes using coevolution and STRING"),
                                   h4("It is recommended to begin with the user guide in the 'Help' tab in case you are new to using this application.",style="text-align:center;"),
                                   br(),br(),
                                   h2("Let's start!"),
                                   actionButton("toggleSidebar", label = "Show search sidebar",
                                                style="background-color:#e6bf70; font-size:22px;")
                               ) 
                           )), #End of fluidRow
                         hr(),
                         h6(HTML("All rights reserved to the ","<a href='https://www.tabach-lab.com/' target='_blank'>Tabach lab.</a>"),style="text-align:center;"),
                         h6("Developed by Christina Canavati of the Tabach lab. Maintained by Christina Canavati", style="text-align:center;")
                         
                ), #End of tabPanel
                
                
                tabPanel(title ="Step 1: Semantic Similarity-based Prioritization", id="SS_tab", value = "Step 1: Semantic Similarity-based Prioritization", 
                         
                         conditionalPanel(
                           condition = "input.submitbutton != 0", 
                           
                           fluidRow(
                             br(),
                             
                             column(12,
                                    
                                    h3("No Candidate gene identified? Click on the button below to navigate to coevolution and STRING based prioritization"),
                                    br(),
                                    actionButton("indirect_button",'Click Here to navigate to coevolution and STRING based prioritization' ,
                                                 icon("dna"), 
                                                 style="color: #fff; background-color: #2d67cc; border-color: #2e6da4; font-size: 22px")
                                    , align = "center"
                                    , style = "margin-bottom: 20px;"
                                    , style = "margin-top: -10px;"
                             ),
                             
                             br(),br(),
                             
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             ),
                             
                             br(),br(),
                             
                             box(width = 12,
                                 title=span("Semantic Similarity-based Results:", style = "color: white;"),
                                 status = "warning",
                                 solidHeader = TRUE, 
                                 h4(HTML("<b>Click on the desired hyperlinked gene to view best matching phenotypes</b>"),style="color:orange;"), 
                                 h5(HTML("<b>Note: if a gene has no phenotype matching your input query, no plot will be displayed</b>"),style="color:orange;"),
                                 br(),
                                 dataTableOutput("tabledata_HPO")), 
                             verbatimTextOutput('contents2'),
                             
                             column(12,
                                    br(), br(),
                                    actionButton("indirect_button2",'Click Here to navigate to coevolution and STRING based prioritization' ,
                                                 icon("dna"), 
                                                 style="color: #fff; background-color: #2d67cc; border-color: #2e6da4; font-size: 22px")
                                    , align = "center"
                                    , style = "margin-bottom: 20px;"
                                    , style = "margin-top: -10px;"
                                    ,br(), br()),
                           )) #End of Conditional Panel
                         
                ), #End of tabPanel
                
                tabPanel(title ="Step 2: Coevolution and STRING-based Prioritization", id="priotization_table_tab", value = "Step 2: Coevolution and STRING-based Prioritization", 
                         
                         conditionalPanel(
                           condition = "input.submitbutton != 0", 
                           
                           br(),br(),
                           
                           fluidRow(
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             ),
                             
                             box(title = span("Barplot", style = "color:white;"),
                                 width = 6,  style = "background-color:#ffffff;",
                                 height = "100%",
                                 status = "primary", 
                                 solidHeader = TRUE, 
                                 collapsible = TRUE, 
                                 plotlyOutput("barplot_box", height = 'auto')
                             ), #End of box  
                             
                             box(width = 6,
                                 title = span("Coevolution and STRING-based gene prioritization results:", style = "color: white;"),
                                 status = "primary", 
                                 solidHeader = TRUE, style ="font-size:16px",
                                 
                                 h4("- The more significant the EvORanker p-value is, the higher number of coevolving and/or STRING interacting genes that are related to the input HPO terms"),
                                 h4(HTML("<b>- Click on the hyperlinked gene name to explore the evolutionary and STRING links between your gene of interest and the HPO-related genes</b>"),style="color:orange;"), 
                                 br(),
                                 dataTableOutput("tabledata")
                             )#End of box
                           )#End of fluidRow
                         )
                ), #End of tabPanel
                
                
                
                
                tabPanel(title = "Per Gene Analyses", id="plots_tab", value = "Per Gene Analyses",
                         
                         conditionalPanel(
                           condition = "input.submitbutton != 0", 
                           
                           fluidRow(
                             box(width = 12, solidHeader = TRUE, collapsible = TRUE,
                                 title = span("Co-evolution and STRING Subetwork", style = "color: white;"),
                                 status = "primary",
                                 h4("Subnetwork showing the top HPO-related genes (above the 85th percentile) within the top 50 genes co-evolving and interacting (STRING) with the query gene."), 
                                 h4("Blue edges indicate co-evolution, red edges indicate STRING interactions"),
                                 h4("The width of the edge correlates with the number of clades the gene is found to co-evolve with the query gene"),
                                 h4("Click the desired nodes in the network to check the HPO terms associated with the desired gene"),
                                 h4(HTML("<b>Hit the 'Download network' button bellow to save network</b>"), style="color:orange;"),
                                 br(),
                                 downloadButton('downloadNetwork', 'Download network',class = "butt"),
                                 tags$head(tags$style(".butt{background-color:#add8e6;} .butt{color: #337ab7;}")),
                                 visNetworkOutput("mynetworkid", height = "1000px", width = "100%")
                                 
                             )   
                           ),
                           
                           fluidRow(
                             box(width = 12, solidHeader = TRUE, collapsible = TRUE,
                                 title = span("Enrichment Analysis", style = "color: white;"),
                                 status = "primary",
                                 h4("Click on the button below to analyze the enriched pathways/ ontologies of all the coevolving and STRING interacting disease genes with your gene of interest with Enrichr!"), 
                                 br(),
                                 actionButton("enrichr_button", span("Analyse in ", span(style="letter-spacing: 0.1em;","En",span("rich", .noWS = c("after", "before"), style="color:red;"),"r")),
                                              style="font-size:14px; border:1px solid #999; border-radius: 2px; position: relative;"),
                                 br(),br(),
                                 h4("Click on the download button below to download the full table of the scored HPO genes within the top 50 coevolving and STRING interacting genes with your gene of interest"), 
                                 br(),
                                 downloadButton("downloadData", "Download Full Data Per Gene",class = "butt"),
                                 tags$head(tags$style(".butt{background-color:#add8e6;} .butt{color: #337ab7;}"))
                                 
                             )
                           )
                           
                         )
                ), #End of tabPanel
                
                
                tabPanel("About",
                         fluidRow(column(10,offset = 1,
                                         HTML( 
                                           '<p><u><span style="font-size: 12.0pt; line-height: 107%;">About:</span></u></p>
                                                            <p>EvORanker was developed by the Prof. Yuval Tabach group for co-evolution-based gene prioritization.</p>
                                                            <p>Our focus is to use comparative genomics and other omics from the STRING database to establish a possible \'missing\' link between a patient\'s candidate gene and a patient\'s phenotypes </p>
                                                            <p>The following are needed as input for EvORanker:</p>
                                                            <ol>
                                                            <li>Patient clinical information as HPO (Human Phenotype Ontology) terms. Please refer to the HPO database to retrieve the HPO terms describing the patient: https://hpo.jax.org/app/ </li>
                                                            <li>Patient candidate genes from your sequencing experiment. Please use HGNC symbols</li>
                                                            </ol>
                                                           
                                                            <p style="tab-stops: 56.2pt;"><strong>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </strong></p>
                                                            <p><strong><span style="font-size: 12.0pt; line-height: 107%;">Cite EvORanker:</span></strong></p>
                                                            <p><strong>Citation will be available once accepted </strong></p>
                                                            <p style="tab-stops: 56.2pt;"><strong>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </strong></p>
                                                            <p><u><span style="font-size: 12.0pt; line-height: 107%;">Additional Reference: </span></u></p>
                                                            <ol>
                                                            
                                                            <li>Tabach, Y., et al., Identification of small RNA pathway genes using patterns of phylogenetic conservation and divergence. Nature, 2013. 493(7434): p. 694-698.</li>
                                                            <li style="margin-bottom: .0001pt;">Tabach Y., et la. Human disease locus discovery and mapping to molecular pathways through phylogenetic profiling. Mol Syst Biol, 2013,9: p. 692.</li>
                                                            <li style="margin-bottom: .0001pt;">Sherill-Rofe, D., et al., Novel players in the homologues recombination repair pathway revealed through local and global co-evolution analysis across hundreds of genomes Accepted, 2018.</li>
                                                            <li style="margin-bottom: .0001pt;">Szklarczyk et al. Nucleic acids research 47.D1 (2018): D607-D613.2</li>
                                                            <li style="margin-bottomL .0001pt:">Greene D, Richardson S, Turro E (2017). “ontologyX: a suite of R packages for working with ontological data.” Bioinformatics, btw763.</li>
                                                            <li style="margin-bottomL .0001pt:">Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Maayan A. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14).</li>
                                                            <li style="margin-bottomL .0001pt:">Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Maayan A.Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90.</li>
                                                            <li style="margin-bottom: .0001pt;">Tsaban Tomer, Doron Stupp, Dana Sherill-Rofe, Idit Bloch, Elad Sharon, Ora Schueler-Furman, Reuven Wiener, and Yuval Tabach. "CladeOScope: functional interactions through the prism of clade-wise co-evolution." NAR Genomics and Bioinformatics 3, no. 2 (2021): lqab024.</li>
                                                            </ol>
                                                            <p>&nbsp;</p>
                                                            <p><u><span style="font-size: 12.0pt; line-height: 107%;">Contact:</span></u></p>
                                                            <p>canavatichristina@gmail.com</p>
                                                            <p><strong>&nbsp;</strong></p>'
                                         )
                         ))        
                ), #End of tabPanel
                
                
                tabPanel("Help",
                         fluidPage(
                           fluidRow(column(10, offset = 1,
                                           tags$iframe(style="height:600px; width:100%; scrolling=yes", src="Help_File.pdf")
                                           #withSpinner(slickROutput("slickr", width="1280px",height = "720px"))
                           )#close column
                           )#close fluidrow
                         ))##closing tabPanel (1) 
                
    )#tabsetPanel
  )#dashboard body
)#dashboard page

####################################
# Server                           #
####################################

server <- function(input, output, session) {
  
  plist <- list()
  glist <- list()
  clist <- list()
  hpolist <- list()
  # Sidebar 
  observeEvent(input$toggleSidebar, {
    if(input$sidebarCollapsed)
      shinyjs::toggleClass(selector = "body", class = "sidebar-collapse")
  })
  
  observeEvent(input$submitbutton,  {
    
    req(input$txt1 != "")
    req(input$txt2 != "")
    #  ifelse(input$txt1 != "" && input$txt2 != "",
    updateTabsetPanel(session = session, "tabs", selected = "Step 1: Semantic Similarity-based Prioritization")
    
  })
  
  observeEvent(input$indirect_button,  {
    
    req(input$txt1 != "")
    req(input$txt2 != "")
    req(input$submitbutton>0)
    #  ifelse(input$txt1 != "" && input$txt2 != "",
    updateTabsetPanel(session = session, "tabs", selected = "Step 2: Coevolution and STRING-based Prioritization")
    #      showModal(modalDialog(size = "m", div("Please provide input", style="font-size:160%")))) 
  })
  
  observeEvent(input$indirect_button2,  {
    
    req(input$txt1 != "")
    req(input$txt2 != "")
    req(input$submitbutton>0)
    #  ifelse(input$txt1 != "" && input$txt2 != "",
    updateTabsetPanel(session = session, "tabs", selected = "Step 2: Coevolution and STRING-based Prioritization")
    #      showModal(modalDialog(size = "m", div("Please provide input", style="font-size:160%")))) 
  })
  # Input Data
  datasetInput <- reactive({  
    
    req(input$txt1 != "")
    req(input$txt2 != "")
    
    showModal(modalDialog(size = "m", div("Computing results, please wait...", style="font-size:160%"), footer=NULL))
    
    # Take HP terms from user and Split query set
    query_set <- unlist(strsplit(str_replace_all(input$txt1, " ", ""),split=','))
    
    valid_check_hpos <- query_set %in% hpo$id
    
    if(length(which(valid_check_hpos=="FALSE"))!=0){
      showModal(modalDialog(size = "m", div(paste("Invalid query, the following input HPO term(s):", paste(query_set[!valid_check_hpos], collapse=", "), "could not be processed, please try again"), style="font-size:160%")))
    }
    
    query_set <- unique(query_set)
    #Take Gene input from user and split
    
    genes <- unlist(strsplit(str_replace_all(input$txt2, "\n", ","),split=','))
    
    comma_vec<-gsub(';|,|_|"|%|$|#|*|@'," ", genes)
    validation_gene_vector<-unlist(strsplit(comma_vec,split=" ",fixed=F))
    validation_gene_vector<-toupper(validation_gene_vector)
    validation_gene_vector<-validation_gene_vector[validation_gene_vector>0]
    
    valid_check_genes <- validation_gene_vector %in% All_genes
    
    if(length(which(valid_check_genes=="FALSE"))!=0){
      showModal(modalDialog(size = "m", div(paste("Invalid query, the following input gene(s):", paste(genes[!valid_check_genes], collapse=", "), "could not be processed, please try again"), style="font-size:160%")))
    }
    
    req(length(which(valid_check_genes=="FALSE"))==0)
    
    genes<-unique(validation_gene_vector)
    
    # Get similarity grid
    sim_grid <- get_sim_grid(ontology=hpo, term_sets=hpo_genes, term_sets2=list(query_set), combine = "average")
    
    # Convert matrix into a dataframe                 
    df.hpo <- as.data.frame(sim_grid)
    df.hpo <- tibble::rownames_to_column(df.hpo, "Gene")
    colnames(df.hpo)[2] <- "Score"
    df.hpo_all <- df.hpo %>% arrange(desc(df.hpo$Score))
    colnames(df.hpo_all)[1] <- "Gene2"
    #df.hpo_all$Score <- round(df.hpo_all$Score,2)
    df.hpo_85 <- df.hpo_all[ df.hpo_all$Score > quantile(df.hpo_all$Score, 0.85 ) , ]
    
    # Generate phenotype similarity score 
    table_HPO <- df.hpo_all %>% filter(Gene2 %in% genes)
    if(nrow(table_HPO) != 0) {
      genes_HPO <- table_HPO$Gene2
      select_hpo_genes <- hpo_genes[genes_HPO]
      sim <- lapply(select_hpo_genes, function(gene_terms) get_term_sim_mat(
        ontology=hpo,
        information_content=IC,
        row_terms=gene_terms,
        col_terms=query_set))
      sim_DT <- lapply(sim, function(x) {data.table(x, keep.rownames = TRUE)})
      sim_DT_Total <- bind_rows(sim_DT, .id = "Gene")
      gathered_DT <- gather(sim_DT_Total, query, value, -c(rn,Gene)) 
      gathered_DT <- gathered_DT %>% group_by(query,Gene) %>%  summarize(value = max(value))  %>% left_join(gathered_DT,by=c("value","query","Gene"))
      gathered_DT <- data.table(gathered_DT)[value != 0]
      named_DT <- gathered_DT %>% mutate(HPO_Name_query = hpo$name[query], HPO_Name_gene = hpo$name[rn]) 
      named_DT$HPO_Name_gene <- paste0(named_DT$HPO_Name_gene, " (", named_DT$rn, ")")
      named_DT$HPO_Name_query <- paste0(named_DT$HPO_Name_query, " (", named_DT$query, ")")
      named_DT <- data.table(named_DT)[value > 0.3]
    } else {
      named_DT <- data.frame(query=0,Gene=0,value=0,rn=0,HPO_Name_query=0, HPO_Name_gene=0)
      named_DT <- named_DT[NULL,]
      
    }
    
    ###################################3      
    
    savethis <- data.frame(Gene = 0, EvORanker_p_value = 0)
    savethis <- savethis[NULL, ]
    
    for(gene in genes){
      
      print(gene)
      
      #Best matching phenotypes plot
      if(nrow(data.table(named_DT)[Gene==gene]) != 0){
        
        hpolist[[gene]] <- data.table(named_DT)[Gene==gene] %>% ggplot(aes(HPO_Name_query, HPO_Name_gene, fill = value)) + 
          geom_tile()  + theme_base() +
          theme( plot.title = element_text(face = "bold", size = 16), axis.text=element_text(size=12),
                 axis.text.x = element_text(angle = 25 , hjust = 1)) +  
          scale_fill_gradientn(limits  = range(0, 1),
                               colours =  myPal2(30),
                               values = scales::rescale( c(0,1), to = c(0, 1), from = c(0, 1))) + xlab("Query HPO Terms") +
          ylab("Gene-associated HPO terms") + ggtitle(paste("Best matching phenotypes for", gene)) 
      }
      
      print("Calculating EvORanker Score")
      
      # Defining all clades, retrieving the HPO genes: 
      Gene_Rank_File <- All_Clades_File[Gene1 == gene]
      Gene_Rank_File <- left_join(Gene_Rank_File, ind, by ="Clade_Index")
      Gene_Rank_File <- left_join(Gene_Rank_File, df.hpo_all) %>% mutate(group = "Coevolution")
      Coevolved_genes <- Gene_Rank_File[,-1]
      colnames(Coevolved_genes) <- c("Gene1", "Gene2", "Coevolutionary Rank", "Clade", "Semantic Similarity Score", "Source")
      Gene_Rank_File <- Gene_Rank_File[,-c(1,4,5)]
      
      # Integrate STRING
      string_gene <- string[P1 == gene]
      colnames(string_gene)[3] <- "Gene2"
      string_gene <- inner_join(string_gene, df.hpo_all)
      string_gene <- data.table(string_gene)
      
      STRING_genes <- string_gene %>% mutate(group = "STRING")
      colnames(STRING_genes) <- c("STRING Combined Score", "Gene1", "Gene2", "Semantic Similarity Score", "Source")
      Combined_DT <- bind_rows(Coevolved_genes, STRING_genes)
      Combined_DT[is.na(Combined_DT)] <- "-"
      clist[[gene]] <- Combined_DT
      
      string_gene <- string_gene[,-1]
      colnames(string_gene) <- c("Gene1", "Gene2", "Score")
      if(nrow(string_gene) > 0){
        string_gene <- string_gene %>% mutate(group = "STRING")
        Gene_Rank_File <- rbind(Gene_Rank_File, string_gene)}
      glist[[gene]] <- Gene_Rank_File %>% filter(Gene2 %in% df.hpo_85$Gene2)
      
      dt <- data.frame(Gene = 0, group = 0, p_value =0)
      dt <- dt[NULL,]
      for(i in unique(Gene_Rank_File$group)){
        
        Group <- Gene_Rank_File[group == i]
        if(nrow(Group) == 0){
          next
        }
        test <- ks.test(Group$Score, df.hpo_all$Score, alternative = "less", exact = FALSE)
        p_value = test$p.value
        tmp <- data.frame(Gene = gene, group = i, p_value = p_value)
        
        dt <- rbind(dt, tmp)
      }
      
      
      #combined_p_value <- min(savethis$p_value)
      combined_p_value = combine.test(dt$p_value, method = "fisher")
      if(length(combined_p_value) == 0){combined_p_value = 1}
      #combined_p_value = pchisq((sum(log(savethis$p_value))*-2), df=length(savethis$p_value)*2, lower.tail=F)
      #combined_p_value <- mean(savethis$p_value)
      dt_tmp <- data.frame(Gene = gene, EvORanker_p_value = combined_p_value)
      #dt_tmp <- data.frame(Gene = gene, p_value = combined_p_value)
      
      savethis <- rbind(savethis, dt_tmp)
      
      #Barplot:
      
      output$barplot_box <- renderPlotly({ 
        p <-  ggplot(savethis, aes(-log10(EvORanker_p_value), reorder(Gene,-log10(EvORanker_p_value)),  fill = Gene)) + #, alpha = 0.8
          geom_bar(stat = "identity", width = 0.5) + #color = "black"
          theme_minimal() +
          theme( axis.text.y = element_text(size = 8), legend.position = "none") + xlab("-log10 EvORanker p-value") + ylab("Genes") 
        ggplotly(p, tooltip="text", height = ifelse(length(genes) %in% c(1:20), 500, 1200))
      })
      
    }
    
    
    removeModal()
    
    savethis[is.na(savethis)] <- 0    
    savethis <- savethis %>% arrange(EvORanker_p_value) %>% mutate(EvORanker_p_value = signif(EvORanker_p_value, digits = 3)) 
    #savethis <- savethis %>% mutate(p_value = ifelse(EvORanker_Score == 0, NA, p_value))
    names(hpolist) <- gsub('-','_',names(hpolist))
    names(clist) <- gsub('-','_',names(clist))
    names(glist) <- gsub('-','_',names(glist))
    result <- list(table_HPO, hpolist, savethis, plist, glist, clist)
    return(result)
    
    
  })
  
  #Outout semantic similarity table
  output$tabledata_HPO <-  DT::renderDataTable({
    req(input$submitbutton)
    if(input$submitbutton>0){
      isolate(datasetInput()[[1]]) 
    } 
    tab1 <- isolate(datasetInput()[[1]])
    shinyInput <- function(FUN, len, id, label, ...) {
      inputs <- character(len)
      
      for (i in seq_len(len)) {
        label <- tab1$Gene2[i]
        inputs[i] <- as.character(FUN(paste0(label,"_", i),label=label, ...))
      }
      inputs
    }
    tab1$Score <- round(tab1$Score,2)
    tab1 = tab1 %>% mutate(Gene2 = shinyInput(actionLink, nrow(tab1),  label, label = Gene2, onclick = 'Shiny.setInputValue(\"select_button1\", this.id, {priority: \"event\"})'))
    tab1 <- datatable(tab1,  callback = JS("
var tips = ['Gene Names', 'Gene Names', 'The higher the Simantic similarity score, the closer the gene is to the phenotype'],
    header = table.columns().header();
for (var i = 0; i < tips.length; i++) {
  $(header[i]).attr('title', tips[i]);
}
"), escape = FALSE 
                      , filter="top", selection="multiple"
                      , colnames = display_colnames_HPO
                      , extensions = c("Buttons")
                      , options = list( 
                        
                        columnDefs = list(list(
                          className = 'dt-center', targets = "_all")),
                        
                        dom = "Blfrtip"
                        , buttons = 
                          list("copy", list(
                            extend = "collection"
                            , buttons = c("csv", "excel", "pdf")
                            , text = "Download"
                          ) )))
    
    tab1
  }) 
  
  
  observeEvent(input$select_button1, {
    selectedRow <- as.numeric(strsplit(input$select_button1, "_")[[1]][2])
    gene_name <- str_remove(input$select_button1, "_.+")
    check_names <- isolate(datasetInput()[[2]])
    gene_name <- gsub('-','_',gene_name)
    if(gene_name %in% names(check_names)){
      showModal(modalDialog(plotOutput("plot_HPO"), size = "l"))
      output$plot_HPO <- renderPlot(isolate({eval(parse(text = paste0("check_names$", gene_name)))
      })
      )
    } else {
      showModal(modalDialog(size = "m", div("No best matching phenotypes to display"), style="font-size:160%"))
      
      
    } 
    
  })
  
  #Output co-evolution based prioritzation table: 
  output$tabledata <-  DT::renderDataTable({
    #req(input$indirect_button)
    if(input$submitbutton>0){
      isolate(datasetInput()[[3]]) 
    } 
    tab <- isolate(datasetInput()[[3]])
    shinyInput <- function(FUN, len, id, label, ...) {
      inputs <- character(len)
      
      for (i in seq_len(len)) {
        label <- tab$Gene[i]
        inputs[i] <- as.character(FUN(paste0(label,"_", i),label=label, ...))
      }
      inputs
    }
    
    tab = tab %>% mutate(Gene = shinyInput(actionLink, nrow(tab),  label, label = Gene, onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})'))
    tab <- datatable(tab,  callback = JS("
var tips = ['Gene Names', 'Gene Names', 'The more significant the EvORanker p-value, the higher the number of coevolving and STRING interacting HPO-related genes'],
    header = table.columns().header();
for (var i = 0; i < tips.length; i++) {
  $(header[i]).attr('title', tips[i]);
}
"), escape = FALSE 
                     , filter="top", selection="multiple"
                     , colnames = display_colnames
                     , extensions = c("Buttons")
                     , options = list( 
                       "pageLength" = 25,
                       columnDefs = list(list(
                         className = 'dt-center', targets = "_all")),
                       
                       dom = "Blfrtip"
                       , buttons = 
                         list("copy", list(
                           extend = "collection"
                           , buttons = c("csv", "excel", "pdf")
                           , text = "Download"
                         ) )))
    
    tab
  }) 
  
  
  #Generate coevolution network
  
  observeEvent(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    updateTabsetPanel(session = session, "tabs", selected = "Per Gene Analyses")
    gene_name <- str_remove(input$select_button, "_.+")
    showModal(modalDialog(size = "m", div("Constructing co-evolution network...", style="font-size:160%"), footer=NULL))
    
    gene_name = gsub('-','_',gene_name)
    rank50 <- isolate(eval(parse(text = paste0("datasetInput()[[5]]$", gene_name))))
    #   if(nrow(rank50) >= 100){rank50 <- rank50[Rank <= 20]}
    net <- All_Clades_File[Gene1 %in% c(gene_name, rank50$Gene2) & Gene2 %in% c(gene_name,rank50$Gene2)]
    if(nrow(net) > 100){net <- net[Rank <= 20]}
    edges <- data.frame(from = net$Gene1, to = net$Gene2, Clade = net$Clade)
    edges <- edges[!duplicated(apply(edges,1,function(x) paste(sort(x),collapse=''))),]
    edges <- edges %>% group_by(from, to) %>% summarize(width = n())
    edges <- edges %>% mutate(group = "Coevolution")
    str_net <- string[P1 %in% c(gene_name, rank50$Gene2) & P2 %in% c(gene_name,rank50$Gene2)]
    str_net <- str_net %>% mutate(Rank = frank(-combined_score, ties.method = "dense"))
    if(nrow(str_net) > 100){str_net <- str_net[Rank <=20]}
    str_net <- str_net %>% select(P1, P2) %>% mutate(group = "STRING")
    nodes <- data.frame(id = unique(c(net$Gene1,net$Gene2, str_net$P1, str_net$P2)))
    nodes <- nodes %>% mutate(label = id)
    colnames(str_net) = c("from", "to", "group")
    edges <- bind_rows(edges, str_net)
    edges[is.na(edges)] <- 4
    edges <- edges %>% mutate(color = ifelse(group == "Coevolution", "blue", "red"))
    ledges <- data.frame(color = c("blue", "red"), 
                         label = c("Coevolution", "STRING"), font.align = "top" )
    # lnodes <- data.frame(color = c("darkblue", "orange"), 
    #                       label = c("HPO-related Gene", "Query Gene"))
    edges <- edges[!duplicated(apply(edges,1,function(x) paste(sort(x),collapse=''))),]
    nodes <- nodes %>% mutate(group = ifelse(id == gene_name,"Query Gene", "HPO-related Gene"))
    nodes <- nodes %>% mutate(font.size = 18)#, color = ifelse(group == "Query Gene", "orange", "darkblue"))
    nodes <- nodes %>% mutate(Gene_symbol = label)
    nodes <-  left_join(nodes, df_hpo[,c("Gene_ID", "Gene_symbol")], by = "Gene_symbol") %>% unique()
    nodes <- nodes %>% mutate(title = ifelse(group!="Query Gene",paste0('<a target="_try" href = "https://hpo.jax.org/app/browse/gene/',nodes$Gene_ID,'\">Check Associated HPO terms</a>'), "Query Gene"))
    final_network <- visNetwork(nodes, edges) %>% visLegend(addEdges = ledges)#,addNodes = lnodes) 
    
    removeModal()
    
    output$mynetworkid <-  renderVisNetwork({final_network %>% visOptions(highlightNearest = TRUE) %>%
        visPhysics(enabled = TRUE) %>% visEdges(smooth = TRUE)})
    
    
    # get position info
    # format positions
    nodes_positions <- reactive({
      positions <- input$network_positions
      if(!is.null(positions)){
        nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
        nodes_positions$id <- names(positions)
        nodes_positions
      } else {
        NULL
      }
    })
    
    
    output$downloadNetwork <- downloadHandler(
      filename = function() {
        paste('network-', Sys.Date(), '.html', sep='')
      },
      content = function(con) {
        nodes_positions <- nodes_positions()
        if(!is.null(nodes_positions)){
          nodes_save <- merge(nodes, nodes_positions, by = "id", all = T)
        } else  {
          nodes_save <- nodes
        }
        
        visNetwork(nodes = nodes_save, edges = edges, height = "800px") %>%
          visOptions(highlightNearest = TRUE) %>% visExport() %>%
          visPhysics(enabled = TRUE) %>% visEdges(smooth = TRUE) %>% visSave(con)
      }
    )
    
    #ENRICHR
    observeEvent(input$enrichr_button, {
      req(input$select_button)
      
      C_Genes <- isolate(unique(c(eval(parse(text = paste0("datasetInput()[[6]]$", gsub('-','_',gene_name))))$Gene2), gene_name))
      genes_str = paste0(C_Genes, collapse = "\n")
      description_str = paste0("Coevolved and STRING-interacting genes with ", gene_name)
      r <- POST(ENRICHR_ADDLIST_URL,
                body = list("list" = genes_str, "description" = description_str)
      )
      
      if(status_code(r)==200) {
        response <- fromJSON(httr::content(r, as = "text"))  # Collect response
        enrichr_permalink <- paste0(ENRICHR_DATASET_URL,
                                    response$shortId[1])
        js$browseURL(enrichr_permalink)
        Sys.sleep(1)
      }
      
      
    })
    gene_name = gsub('-','_',gene_name)
    content = eval(parse(text = paste0("datasetInput()[[6]]$", gene_name)))
    output$downloadData <- downloadHandler(
      filename = paste0("Top-50-Coevolved-and-STRING-Genes-with-",gene_name,".csv"),
      content = function(file) {
        write.csv(content, file, row.names = FALSE)
      }
    )
  })
  
  
  
  
}
####################################
# Create the shiny app             #
####################################
shinyApp(ui = ui, server = server)
