##### Install the packages


source("http://bioconductor.org/biocLite.R")  # Installer for bioconductor packages
biocLite("org.Hs.eg.db")   # Database for converting Entrez Gene ID to NAMEs/SYMBOLs
biocLite("DOSE")  # DOSE package for Disease Ontology
install.packages("rentrez")  # Package used to query Pubmed using EUtils to get information
install.packages("data.table")
install.packages("RISmed") # Package similar to rentrez for querying PubMed
install.packages("pubmed.mineR") # Pubtator Miner to extract chemicals from abstract text

##### Loading the packages


library(DO.db)   # Loading the disease ontology package
library(rentrez)   # Loading the rentrez package
library(data.table)
library(org.Hs.eg.db) # To load the database from bioconductor for genes
library(RISmed) # For EUtilsSummary function for querying pubmed
library(pubmed.mineR)

##### STEP 1 : Go through the whole ontology and get all the DOIDs that belong to cancer


# Recursive function to develop a table based on the ancestor which in this case is DOID 162: cancer
rec_getTable <- function(cur_doid, this_name_preceding, this_preceding_doid, this_level) {
  x <- get (cur_doid,xx)   # get the object based on DOID
  disease_name=Term(x)   # getting the disease name from the object
  doid=DOID(x)   # getting the DOID based on the object
  level = this_level     # getting the levels below the ancestor (which is level 0)
  preceding_doid= this_preceding_doid   # getting the DOID of the preceding disease in the ontology
  name_preceding = this_name_preceding   # getting the disease name of the preceding disease in the ontology
  disease_syn=Term(x)    # getting the synonyms of the disease
  # Combining the data with the data table to populate the fields
  dt <<- rbindlist(list(dt, list(disease_syn, disease_name, doid, level, preceding_doid, name_preceding)))  
  # for each children of the current disease, run this recursive function to populate the data table
  x_child <- get(cur_doid, x_children)
  for (i in 1:length(x_child)) {
    if (!(is.na(x_child[i]))) {
      rec_getTable(x_child[i], disease_name, cur_doid, this_level + 1)
    }
  }
}
# Getting the Disease Ontology IDs (DOIDs) of all the cancer terms (starting from DOID 162 "cancer")
xx <- as.list(DOTERM) # getting the list of disease terms from the Ontology
x_children <- as.list(DOCHILDREN) # getting the list of children of disease
# Data table that contains the resulting information
dt <- data.table(disease_syn=vector(), disease_name=vector(), doid=vector(), level=vector(), preceding_doid=vector(), name_preceding = vector())
rec_getTable("DOID:162","","",0) # Using function recTable to develop the table used to store the cancer disease information


# (TO COMMENT) when implementing on all
doid_in_question <- "DOID:3908"

##### STEP 2 : For each DOID get the Gene IDs that correspond to that

setkey(dt,doid)
doid_keys <- dt$doid  

doid2gene_list <-  vector(mode="list", length=length(doid_keys))
names(doid2gene_list) <- doid_keys
for (i in 1:length(doid_keys))  {
  search_term <- paste(dt[doid_keys[i]]$disease_name, "AND Humans[MeSH Terms]")
  if (doid_keys[i]==doid_in_question) {   # (TO COMMENT) when implementing on all
    print(search_term)
    another_r_search <- entrez_search(db="gene", term=search_term, retmax=10000)  
    doid2gene_list[[doid_keys[i]]] <- another_r_search$ids
  }
  
}



##### STEP 3 : For each of the gene ID get the PMIDs that correspond to that gene
# This can do through using org.Hs.egPMID


# First get the unique gene list from doid2gene_list
v <- c()
for (j in 1:length(doid2gene_list)) {   
  for (i in 1:NROW(doid2gene_list[[j]])) {
    if (!(length(doid2gene_list[[j]])== 0))
      v <- c(v, doid2gene_list[[j]][i])
  }
}
unique_gene_list <- unique(v)

# Then create a list of lists that stores the PMID for each gene
# Getting a list of PMID mapped to entrez gene IDs (how about using eUtil instead?)
x <- org.Hs.egPMID
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
gene2PMID_list <-  vector(mode="list", length=length(unique_gene_list))
names(gene2PMID_list) <- unique_gene_list

for (i in 1:length(unique_gene_list))  {
  try(gene2PMID_list[[unique_gene_list[i]]] <- get(unique_gene_list[i], xx),silent=T)
}




##### STEP 4 :  For each DOID get the PMIDs that correspond to that disease
# entrez_search db = pubmed search by disease
# Choose one disease (NSCLC)  (NEXT TIME CONVERT TO LOOP)


setkey(dt,doid)
doid_keys <- dt$doid  
doid2PMID_list <-  vector(mode="list", length=length(doid_keys))
names(doid2PMID_list) <- doid_keys
search_term <- paste(dt[doid_in_question]$disease_name, "AND Humans[MeSH Terms]")
search_records <- EUtilsSummary(search_term, type="esearch",db = "pubmed",retstart=0, retmax=99999)
doid2PMID_list[[doid_in_question]] <- search_records@PMID


# temp to delete: Write to table
doidpmid_dt <- data.table(doid=vector(), pmid=vector())
for (i in 1:length(doid2PMID_list[[doid_in_question]])) {
  pmid <- doid2PMID_list[[doid_in_question]][i]
  doidpmid_dt <- rbindlist(list(doidpmid_dt, list(doid_in_question, pmid)))  
}
write.table(doidpmid_dt,"doidpmid_dt.csv",sep=",")
##### STEP 5 : Find the intersection between the PMIDs that belong to the genes and the PMIDs that belong to the diseases


DGC_list <-  vector(mode="list", length=length(doid_keys))  # disease-gene-connection list 
names(DGC_list) <- doid_keys
DGC_list[[doid_in_question]] <- vector(mode="list",length = length(doid2gene_list[[doid_in_question]]))
names(DGC_list[[doid_in_question]]) <- doid2gene_list[[doid_in_question]]
# From doid2gene_list find all the genes associated with particular disease
for (i in 1:length(doid2gene_list[[doid_in_question]])) {
  # For each gene in doid, run through the gene2PMID_list to get PMIDs for each gene
  gene_pmid_list <- gene2PMID_list[[doid2gene_list[[doid_in_question]][i]]]
  # For the doid get the list of PMID
  disease_pmid_list <- doid2PMID_list[[doid_in_question]]
  # find the intersection of these 2 PMID lists
  DGC_list[[doid_in_question]][[i]] <- intersect(gene_pmid_list, disease_pmid_list) #intersection
  
}

# temp: writing dgc list to data table.
dgc_dt <- data.table(doid=vector(), geneid=vector(), pmid=vector())
namesOfGenes <- names(DGC_list[[doid_in_question]])
for (i in 1:length(namesOfGenes)) {
  if (!(length(DGC_list[[doid_in_question]][[namesOfGenes[i]]])==0)) {
    abc <- cbind(doid=doid_in_question, geneid = namesOfGenes[i], pmid = DGC_list[[doid_in_question]][[namesOfGenes[i]]])
    dgc_dt <- rbindlist(list(dgc_dt, data.frame(abc))) 
  }
}
write.table(dgc_dt, "dgc_dt.csv",sep=",")

    # confidence = no. in intersection / total no. of pmids for doid2pmid_list


##### STEP 6 : For each unique gene that is covered find the other information such as the gene symbol, gene name and perhaps alias
##### this one you need the org.Hs.eg.db
# Create a data table to store the gene information based on ID
gene_dt <- data.table(gene_id=vector(), gene_symbol=vector(), gene_name=vector())


# Getting the list of gene symbols mapped to entrez gene ID
g_sym <- org.Hs.egSYMBOL 
mapped_genesym <- mappedkeys(g_sym)
g_sym_list <- as.list(g_sym[mapped_genesym])
# Getting the list of gene names mapped on entrez gene ID
g_name <- org.Hs.egGENENAME
mapped_genename <- mappedkeys(g_name)
g_name_list <- as.list(g_name[mapped_genename])

# Loop through all the genes in unique_gene_list
for (i in 1:length(unique_gene_list)) {
  try ( {
    
  gene_symbol <- get(unique_gene_list[i], g_sym_list)
  gene_name <- get(unique_gene_list[i], g_name_list)
  gene_dt <- rbindlist(list(gene_dt, list(unique_gene_list[i], gene_symbol, gene_name)))  
  } , silent = T)
}
write.table(gene_dt,"gene_dt.csv",sep=",")   # write this to csv (change in the future to make it automatic)

##### STEP 7 : For each PMID found in the intersect you get the chemicals (using PubMed Miner)
# Get unique list of PMID from intersected list
v <- c()
for (i in 1:length(doid2gene_list[[doid_in_question]])) {
  if (!(length(DGC_list[[doid_in_question]][[i]])== 0))
     v <- c(v, DGC_list[[doid_in_question]][[i]])
}
unique_PMID_list <- unique(v)
PMID2Chemical_list <- vector(mode="list", length=length(unique_PMID_list)) 
names(PMID2Chemical_list) <- unique_PMID_list
for (i in 1:length(unique_PMID_list)) {
  extracted_info <- pubtator_function(unique_PMID_list[i])
  try(PMID2Chemical_list[[unique_PMID_list[i]]] <- extracted_info$Chemicals, silent=T)
}


PMID2Chemical_dt <- data.table(pmid=vector(), chemical=vector())

# Prepare to make it into data table
for (i in 1:length(unique_PMID_list)) {
  if (!(length(PMID2Chemical_list[[unique_PMID_list[i]]])==0)) {
    for (j in 1:length(PMID2Chemical_list[[unique_PMID_list[i]]])) {
      PMID2Chemical_dt <- rbindlist(list(PMID2Chemical_dt, list(unique_PMID_list[i], PMID2Chemical_list[[unique_PMID_list[i]]][j])))         
    }
  }
}
write.table(PMID2Chemical_dt,"PMID2Chemical.csv",sep=",")


# right now have a gene
# then have to find the db related to that
all_the_links <- entrez_link(dbfrom='gene', id=1956, db='all')
entrez_summary(db="gene_pubmed_rif",id=26344197)
##### STEP 8 : Storing the information somewhere to be able to visualize



##### Further Improvements
# 1.  Get the conclusion of each PMID
# 2.  Do document categorization to get the type of document (evidence, etc)
# 3.  Do sentiment analysis to get the treatment success rate
