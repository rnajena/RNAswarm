#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tximport)
library(tximportData)

# Import the data for NA_NP
NA_NP <- read_csv('../results/202108/20210824/NA_NP_interactions.csv')

# Run DEseq2 on it

# Generate circus plots





NA_NP_dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")