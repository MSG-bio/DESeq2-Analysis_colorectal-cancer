# DESeq2 Differential Expression Analysis
# Dataset: GSE142279 (colorectal cancer, 20 tumour vs 20 adjacent normal)
# Author: MSG-bio
# Date: May 2026


#Let's load the library 
library(DESeq2)

#Other important packages

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)
library(GEOquery)
library(org.Hs.eg.db)
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

metadata.mod <- metadata %>%
  select(1,2,36) %>%
  rename(tissue = 'tissue:ch1') 

#Reshape the data

#For analysis of gene long format is best. For DEseq2 we need wider format.
#So gene ID = rows and samples = columns

dat <- as.data.frame(dat)
rownames(dat) <- dat$GeneID
dat$GeneID <- NULL
dat[1:5, 1:5]

head(metadata.mod)
rownames(metadata.mod) <- metadata.mod$geo_accession

#Make sure that all the column names in dat correspond to row names in metadata.mod
#Make sure that they are in the same order
#we can use this to test if they match all(colnames(dat) %in% rownames(metadata.mod))
all(colnames(dat) == rownames(metadata.mod))

#Construct a DEseq data set object

dds <- DESeqDataSetFromMatrix(countData = dat, colData = metadata.mod, design = ~ tissue)

#Pre filter : remove rows with low gene counts
#keeping rows that have at least 10 genes total

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dim(dds)

#Set the factor level, our reference level

dds$tissue <- relevel(dds$tissue, ref = "adjacent normal")

#Note: collapse technical replicates if needed

table(metadata.mod$title) # this data has none
# if it did: dds <- collapseReplicates(dds, groupby = dds$title)

#Run Deseq
dds <- DESeq(dds)

#Results

res <- results(dds)
summary(res)
res.sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
nrow(res.sig)
table(res.sig$log2FoldChange > 0)

#Visualization

# convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys      = rownames(res),
                       column    = "SYMBOL",
                       keytype   = "ENTREZID",
                       multiVals = "first")

#Volcano Plot
png("Volcano_plot_GSE142279_final.png", width = 10, height = 8, units = "in", res = 300)

EnhancedVolcano(res,
                lab      = gene_symbols,
                x        = 'log2FoldChange',
                y        = 'padj',
                pCutoff  = 0.05,
                FCcutoff = 1,
                title    = 'Tumour vs Adjacent Normal',
                subtitle = 'Colorectal Cancer — GSE142279')

dev.off()

#heatmap

# order results by adjusted p-value and take top 30
top30 <- order(res$padj)[1:30]

# extract normalized counts for those genes
mat <- counts(dds, normalized = TRUE)[top30, ]

# log2 transform for visualization
mat <- log2(mat + 1)

# swap Entrez IDs for gene symbols as row names
rownames(mat) <- gene_symbols[rownames(mat)]


# create annotation showing which samples are tumor vs normal
annotation <- data.frame(
  Tissue = metadata.mod$tissue,
  row.names = rownames(metadata.mod)
)

png("Heatmap_GSE142279_final.png", width = 10, height = 12, units = "in", res = 300)

pheatmap(mat,
         annotation_col = annotation,
         show_colnames  = FALSE,
         fontsize_row   = 8,
         main           = "Top 30 DE Genes — GSE142279")

dev.off()