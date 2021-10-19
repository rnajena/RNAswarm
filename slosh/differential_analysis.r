library("DESeq2")

# Import the data for test
input_filepath <- "/home/ru27wav/Projects/gl_iav-splash_freiburg/results/test/wt_mut_interactions.csv"

output_filepath = /home/ru27wav/Projects/gl_iav-splash_freiburg/results/test/deseq2_output.csv

counts <- read.csv(input_filepath, header = FALSE, row.names=1)
conditions <- c("wt", "mut"

# Hack to make the dataframe looks like what it should look to be used as countData DEseq2 input
counts <- counts[,2:4]
colnames(counts) <- samples

# Generate a colData dataframe for DEseq2
col <- data.frame(condition = c("WT","WT","WT"))
row.names(col) <- samples

# Run DEseq2 on it
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col,
                              design= ~ conditions)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Wt NA_NP")

write.csv(as.data.frame(resOrdered), 
          file=output_filepath)
