#!/usr/bin/env Rscript

'Run DESeq2 on two count tables and output a tsv file with the results

Usage:
  run_DESeq2.r --count_table1 <count_table1> --alias1 <alias1> --count_table2 <count_table2> --alias2 <alias2> --output_file <output_file>
  run_DESeq2.r (-h | --help)

Options:
    -h --help                       Show this screen.
    --count_table1 <count_table1>   Path to the count table 1
    --alias1 <alias1>               Alias for the first count table
    --count_table2 <count_table2>   Path to the count table 2
    --alias2 <alias2>               Alias for the second count table
    --output_file <output_file>     Path to the output file

' -> doc

library(docopt)
arguments <- docopt(doc)

# This script runs DESeq2 on two count tables and outputs a tsv file with the results
library("DESeq2")

# Import the data for test, the input is a tsv file with the counts
counts1 <- read.table(arguments$count_table1, header = 1, row.names=1)
counts2 <- read.table(arguments$count_table2, header = 1, row.names=1)


# samples vector is the header of the count table but only the base name of the sample
samples1 <- colnames(counts1)
samples2 <- colnames(counts2)


# Generate a colData dataframe for DEseq2 and merge the two count tables
col <- data.frame(conditions = c(rep("wt", length(samples1)), rep("mutant", length(samples2))))
row.names(col) <- c(samples1, samples2)
counts <- merge(counts1, counts2, by="row.names")
# Use the first column as the row names
rownames(counts) <- counts[,1]
counts <- counts[,-1]


# Run DEseq2 on it and make sure that counts1/sample1 are considered as the reference
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col,
                              design= ~ conditions)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
# Order the results by the pvalue
resOrdered <- res[order(res$pvalue),]


# Write the results to a tsv file, include the rownames column and name it as "interaction_id"
write.table(resOrdered, file=arguments$output_file, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


