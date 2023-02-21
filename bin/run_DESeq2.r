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

# Import the data for test
counts1 <- read.csv(arguments$count_table1, header = 1, row.names=1)
counts2 <- read.csv(arguments$count_table2, header = 1, row.names=1)


# samples vector is the header of the count table but only the base name of the sample
samples1 <- colnames(counts1)
samples2 <- colnames(counts2)


# Generate a colData dataframe for DEseq2 and merge the two count tables
col <- data.frame(conditions = c(rep(arguments$alias1, length(samples1)), rep(arguments$alias2, length(samples2))))
row.names(col) <- c(samples1, samples2)
counts <- merge(counts1, counts2, by="row.names")
counts <- counts[,-1] # remove the first column which is the row names
counts <- as.matrix(counts) # convert to matrix


# Run DEseq2 on it
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col,
                              design= ~ conditions)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
resOrdered <- res[order(res$pvalue),]


# Write the results to a tsv file
write.table(as.data.frame(resOrdered), 
            file=arguments$output_file, 
            sep="\t")
