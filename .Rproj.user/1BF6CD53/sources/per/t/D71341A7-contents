#Install Bioconducter
install.packages("BiocManager")
BiocManager::install("DESeq2")

#Let's load the library 
library(DESeq2)

#Other important packages
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "ggrepel"))

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)
library(GEOquery)
library(org.Hs.eg.db)

#"C:/Users/garci/Desktop/Miriam/Independent learning/R Studies/Portfolio/DEseq2_CC"
library(tidyverse)
library(dplyr)

#read data

dat <- read_tsv(file = 'GSE142279_raw_counts_GRCh38.p13_NCBI.tsv')
dim(dat)
head(dat)
#get metadata

gse <- getGEO(GEO= "GSE142279", GSEMatrix = TRUE)
names(pData(gse[[1]])) 
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
metadata_subset <- select(metadata, c(1,2,10,36))
head(metadata_subset)


metadata.mod <- metadata %>%
  select(1,2,36) %>%
  rename(tissue = 'tissue:ch1') 
  
#Reshape the data

dat.long <- dat %>%
  pivot_longer(cols = starts_with ('GSM'), names_to = 'samples', values_to = 'count')
  

#Join expression data with metadata

dat.joined <- dat.long %>%
  left_join(. , metadata.mod, by = c("samples" = "geo_accession"))

dim(dat)
dim(metadata.mod)

#For analysis of gene long format is best. For DEseq2 we need wider fromet.
#So gene ID = rows and samples = columns


rownames(dat) <- dat$GeneID
dat$GeneID <- NULL
dim(dat)
dat[1:5,1:5]
type(dat)
dat <- as.data.frame(dat)
rownames(dat) <- dat$GeneID
dat$GeneID <- NULL
dat[1:5,1:5]