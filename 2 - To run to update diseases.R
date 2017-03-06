##### Install the packages



##### Loading the packages


library(DO.db)   # Loading the disease ontology package
library(rentrez)   # Loading the rentrez package # Need this
library(data.table)   # Need this
library(org.Hs.eg.db) # To load the database from bioconductor for genes # Need this
library(RISmed) # For EUtilsSummary function for querying pubmed # Need this
library(pubmed.mineR) # Need this
library(RODBC)
library(plyr)



counter = 1
record_number = 1
##### STEP 1 : Set up existing connection to local Access Database
while (counter > 0 ) {
try({
print("Opening Database...")
db <- "PB_Database.accdb"
ch <- odbcConnectAccess2007(db)
doid_table <- sqlQuery(ch, "select * from doid_list")
gene_table <- sqlQuery(ch, "select gene_id from gene_dt")
pmid_table <- sqlQuery(ch, "select distinct pmid from PMID2Chemical")
gene_ids <- gene_table$gene_id
pmid_list <- pmid_table$pmid

sb_data <- subset(doid_table, last_category==1 & is.na(last_update))
counter <- nrow(sb_data)
# CHOOSE the 1st DOID in the table
doid_in_question <- as.character(sb_data$doid[[record_number]])
doid_id_database <- as.character(sb_data$ID1[[record_number]])
disease_name <- as.character(sb_data$disease_name[[record_number]])

##### STEP 2 : For each DOID get the Gene IDs that correspond to that
print("Querying PubMed for genes related to doid")
search_term <- paste(disease_name, "AND Humans[MeSH Terms]")
another_r_search <- entrez_search(db="gene", term=search_term, retmax=9999)  
doid2gene_list <- another_r_search$ids


##### STEP 3 : For each of the gene ID get the PMIDs that correspond to that gene
print("Querying Bioconductor database for PMIds that correspond to gene")
unique_gene_list <- unique(doid2gene_list)
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
print("Querying PubMed for PMIds that correspond to disease")
search_term <- paste(disease_name, "AND Humans[MeSH Terms]")
search_records <- EUtilsSummary(search_term, type="esearch",db = "pubmed",retstart=0, retmax=9999)
doid2PMID_list <- search_records@PMID


##### STEP 5 : Find the intersection between the PMIDs that belong to the genes and the PMIDs that belong to the diseases
DGC_list <- vector(mode="list",length = length(doid2gene_list))
names(DGC_list) <- doid2gene_list
# From doid2gene_list find all the genes associated with particular disease
for (i in 1:length(doid2gene_list)) {
  # For each gene in doid, run through the gene2PMID_list to get PMIDs for each gene
  gene_pmid_list <- gene2PMID_list[[doid2gene_list[i]]]
  # For the doid get the list of PMID
  disease_pmid_list <- doid2PMID_list
  # find the intersection of these 2 PMID lists
  DGC_list[[i]] <- intersect(gene_pmid_list, disease_pmid_list) #intersection
  confidence_level <- length(DGC_list[[i]])/length(doid2PMID_list)
}


##### STEP 6 : For genes that do not already exist in database find the other information such as the gene symbol, gene name and perhaps alias
##### this one you need the org.Hs.eg.db
# Create a data table to store the gene information based on ID
# Find genes that are not in database
new_genes_list <- doid2gene_list[!(doid2gene_list %in% gene_ids)]
gene_dt <- data.table(gene_id=vector(), gene_symbol=vector(), gene_name=vector())

# Getting the list of gene symbols mapped to entrez gene ID
g_sym <- org.Hs.egSYMBOL 
mapped_genesym <- mappedkeys(g_sym)
g_sym_list <- as.list(g_sym[mapped_genesym])
# Getting the list of gene names mapped on entrez gene ID
g_name <- org.Hs.egGENENAME
mapped_genename <- mappedkeys(g_name)
g_name_list <- as.list(g_name[mapped_genename])

# Loop through all the genes in new_genes_list
for (i in 1:length(new_genes_list)) {
  try ( {
  gene_symbol <- get(new_genes_list[i], g_sym_list)
  gene_name <- get(new_genes_list[i], g_name_list)
  gene_dt <- rbindlist(list(gene_dt, list(new_genes_list[i], gene_symbol, gene_name)))  
  } , silent = T)
}


##### STEP 7 : For each PMID found in the intersect that does not already exist in database you get the chemicals (using PubMed Miner)
# Get unique list of PMID from intersected list
print("Querying Pubtator")
v <- c()
for (i in 1:length(doid2gene_list)) {
  if (!(length(DGC_list[[i]])== 0))
     v <- c(v, DGC_list[[i]])
}
unique_PMID_list <- unique(v)
unique_PMID_list <- unique_PMID_list[!(unique_PMID_list %in% as.character(pmid_list))]
PMID2Chemical_list <- vector(mode="list", length=length(unique_PMID_list))

if (length(unique_PMID_list) > 300) {
  print("Skipping current DOID as too many PMID")
  record_number = record_number + 1
  next
}

names(PMID2Chemical_list) <- unique_PMID_list
for (i in 1:length(unique_PMID_list)) {
  extracted_info <- pubtator_function(unique_PMID_list[i])
  try(PMID2Chemical_list[[unique_PMID_list[i]]] <- extracted_info$Chemicals, silent=T)
}


##### STEP 8 : Storing the information to Access database in order to be able to visualize
# FOR UPDATES
# Update doidpmid_dt table in Access Database using 
print("Updating database")
doidpmid_df <- cbind(doid = doid_in_question, pmid = doid2PMID_list)
sqlSave(ch, dat = data.frame(doidpmid_df),tablename = "doidpmid_dt", append = TRUE, rownames = FALSE)

# Update dgc_dt table in Access Database.
dgc_df <- data.frame(doid=vector(), geneid=vector(), pmid=vector())
namesOfGenes <- names(DGC_list)
for (i in 1:length(namesOfGenes)) {
  if (!(length(DGC_list[[namesOfGenes[i]]])==0)) {
    abc <- cbind(doid = doid_in_question, geneid = namesOfGenes[i], pmid = DGC_list[[namesOfGenes[i]]])
    dgc_df <- rbindlist(list(dgc_df, data.frame(abc))) 
  }
}
sqlSave(ch, dat = data.frame(dgc_df),tablename = "dgc_dt", append = TRUE, rownames = FALSE)

# Update gene_dt table in Access DAtabase
sqlSave(ch, dat = data.frame(gene_dt),tablename = "gene_dt", append = TRUE, rownames = FALSE)



# Update PMID2 Chemical List
if (length(unique_PMID_list) > 0 ) {
PMID2Chemical_df <- data.frame(pmid=vector(), chemical=vector())
for (i in 1:length(unique_PMID_list)) {
  if (!(length(PMID2Chemical_list[[unique_PMID_list[i]]])==0)) {
    PMID2Chemical_df <- rbindlist(list(PMID2Chemical_df, data.frame(cbind(pmid = unique_PMID_list[i], chemical = PMID2Chemical_list[[unique_PMID_list[i]]])))) 
  }
}
sqlSave(ch, dat = data.frame(PMID2Chemical_df),tablename = "PMID2Chemical", append = TRUE, rownames = FALSE)
}

# Update Database to indicate that it has been downloaded for particular DOID
updateDOID.df <- sb_data[record_number,]
updateDOID.df$last_update <- Sys.Date()
sqlUpdate(ch, dat = updateDOID.df, tablename = "doid_list",index=c("ID1"))



odbcClose(ch)
}, silent=T)
  
record_number = record_number + 1

}

}

##### Further Improvements
# 1.  Get the conclusion of each PMID
# 2.  Do document categorization to get the type of document (evidence, etc)
# 3.  Do sentiment analysis to get the treatment success rate


