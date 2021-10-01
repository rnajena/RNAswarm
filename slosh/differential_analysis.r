#!/usr/bin/env Rscript
#library("readr")
#library("dplyr")
#library("tximport")
#library("tximportData")
library("DESeq2")

# Import the data for NA_NP
NA_NP_counts <- read.csv('results/202109/20210930/NA_NP_interactions.csv', header = FALSE, row.names=1)
NA_NP_samples <- c("wt01", "wt02", "wt03")
NA_NP_conditions <- c("WT","WT","WT")

# Hack to make the dataframe looks like it should look to be uses as countData DEseq2 input
NA_NP_counts <- NA_NP_counts[,2:4]
colnames(NA_NP_counts) <- NA_NP_samples

# Generate a colData dataframe for DEseq2
NA_NP_col <- data.frame(condition = c("WT","WT","WT"))
row.names(NA_NP_col) <- NA_NP_samples

# Run DEseq2 on it
NA_NP_dds <- DESeqDataSetFromMatrix(countData = NA_NP_counts,
                              colData = NA_NP_col,
                              design= ~ 1)

NA_NP_dds <- DESeq(NA_NP_dds)
resultsNames(NA_NP_dds) # lists the coefficients
NA_NP_res <- results(NA_NP_dds, name="Wt NA_NP")

# Questions I want to answer:
# Which interactions are constant in all samples? And which interactions are not?

# Data visualization 
# Generate circus plots to identify interactions that are differentially presentin our samples


