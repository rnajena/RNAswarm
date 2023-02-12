#!/usr/bin/env Rscript

library("DESeq2")

# Import the data for test
input_filepath <- "/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/8xHA_update0122_CJ_cleaned_counttable.csv"

output_filepath = "/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/8xHA_update0122_CJ_cleaned_counttable_deseq2.csv"

counts <- read.csv(input_filepath, header = 1, row.names=1)
samples <- c("SC35M_WTWT_repli01_0120_trimmed","SC35M_WTWT_repli01_1021_trimmed","SC35M_WTWT_repli02_1120_trimmed","SC35M_WTWT_repli03_1120_trimmed","8-wt_S8_R1_001_trimmed","7wtintakt_S7_R1_001_trimmed","H-4_S14_R1_001_trimmed","SC35M_8xHA_repli01_1021_trimmed","SC35M_8xHA_repli02_1021_trimmed","SC35M_8xHA_repli03_1021_trimmed")

conditions <- c("wt", "8xha")

# Hack to make the dataframe looks like what it should look to be used as countData DEseq2 input
# counts <- counts[,1:6]
colnames(counts) <- samples

# Generate a colData dataframe for DEseq2
col <- data.frame(conditions = c("wt","wt","wt","wt","wt","wt","wt","8xha","8xha","8xha"))
row.names(col) <- samples

# Run DEseq2 on it
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col,
                              design= ~ conditions)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
resOrdered <- res[order(res$pvalue),]

write.csv(as.data.frame(resOrdered), 
          file=output_filepath)
