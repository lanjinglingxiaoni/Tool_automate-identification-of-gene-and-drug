# Installing Packages required

source("http://bioconductor.org/biocLite.R")  # Installer for bioconductor packages
biocLite("org.Hs.eg.db")   # Database for converting Entrez Gene ID to NAMEs/SYMBOLs
biocLite("DOSE")  # DOSE package for Disease Ontology
install.packages("rentrez")  # Package used to query Pubmed using EUtils to get information
install.packages("data.table")  # Data Table for faster processing
install.packages("RISmed") # Package similar to rentrez for querying PubMed
install.packages("pubmed.mineR") # Pubtator Miner to extract chemicals from abstract text
install.packages("RODBC")