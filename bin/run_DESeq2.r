#!/usr/bin/env Rscript

# This script runs DESeq2 on two count tables and outputs a tsv file with the results
library("DESeq2")


# Define command line options
options = list(
    count_table1 = "",
    alias1 = "",
    count_table2 = "",
    alias2 = "",
    output_file = ""
)


# Parse command line options
while(length(commandArgs(TRUE)) > 0) {
    arg = commandArgs(TRUE)[1]
    if(arg == "--count_table1") {
        options$count_table1 = commandArgs(TRUE)[2]
        commandArgs(TRUE) = commandArgs(TRUE)[3:length(commandArgs(TRUE))]
    } else if(arg == "--alias1") {
        options$alias1 = commandArgs(TRUE)[2]
        commandArgs(TRUE) = commandArgs(TRUE)[3:length(commandArgs(TRUE))]
    } else if(arg == "--count_table2") {
        options$count_table2 = commandArgs(TRUE)[2]
        commandArgs(TRUE) = commandArgs(TRUE)[3:length(commandArgs(TRUE))]
    } else if(arg == "--alias2") {
        options$alias2 = commandArgs(TRUE)[2]
        commandArgs(TRUE) = commandArgs(TRUE)[3:length(commandArgs(TRUE))]
    } else if(arg == "--output_file") {
        options$output_file = commandArgs(TRUE)[2]
        commandArgs(TRUE) = commandArgs(TRUE)[3:length(commandArgs(TRUE))]
    } else {
        stop("Unknown argument: ", arg)
    }
}


# Check if required options are provided
if(options$count_table1 == "") {
    stop("Please provide a count table 1 with --count_table1")
}
if(options$alias1 == "") {
    stop("Please provide an alias for count table 1 with --alias1")
}
if(options$count_table2 == "") {
    stop("Please provide a count table 2 with --count_table2")
}
if(options$alias2 == "") {
    stop("Please provide an alias for count table 2 with --alias2")
}
if(options$output_file == "") {
    stop("Please provide an output file with --output_file")
}


# Import the data for test
counts1 <- read.csv(options$count_table1, header = 1, row.names=1)
counts2 <- read.csv(options$count_table2, header = 1, row.names=1)


# samples vector is the header of the count table but only the base name of the sample
samples1 <- colnames(counts1)
samples2 <- colnames(counts2)


# Generate a colData dataframe for DEseq2 and merge the two count tables
col <- data.frame(conditions = c(rep(options$alias1, length(samples1)), rep(options$alias2, length(samples2))))
row.names(col) <- c(samples1, samples2)
counts <- merge(counts1, counts2, by="row.names")
counts <- counts[,-1] # remove the first column which is the row names


# Run DEseq2 on it
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col,
                              design= ~ conditions)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
resOrdered <- res[order(res$pvalue),]

# write.csv(as.data.frame(resOrdered), 
#           file=output_filepath)

# Write the results to a tsv file
write.table(as.data.frame(resOrdered), 
            file=options$output_file, 
            sep="\t")
