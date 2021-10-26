library("DESeq2")

# Import the data for test
input_filepath <- "/home/gabriellovate/Bierinf/Projects/gl_iav-splash_freiburg/results/IAV_wt_vs_mut/wt_mut_interactions.csv"

output_filepath = "/home/gabriellovate/Bierinf/Projects/gl_iav-splash_freiburg/results/IAV_wt_vs_mut/deseq2_output_wt_mut_interactions.csv"

counts <- read.csv(input_filepath, header = 1, row.names=1)
samples <- c("wt01", "wt02", "wt03", "mut01", "mut02", "mut03")
conditions <- c("wt", "mut")

# Hack to make the dataframe looks like what it should look to be used as countData DEseq2 input
counts <- counts[,1:6]
colnames(counts) <- samples

# Generate a colData dataframe for DEseq2
col <- data.frame(conditions = c("wt", "wt","wt","mut","mut","mut"))
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
