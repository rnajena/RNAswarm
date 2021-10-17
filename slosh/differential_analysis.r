library("DESeq2")

# Import the data for NA_NP
args <- commandArgs(trailing = TRUE)
input_filepath <- args[1]

output_filepath = sprintf("%s_deseq2.csv" ,substr(input_filepath, 0, (nchar(input_filepath) - 4)))

counts <- read.csv(input_filepath, header = FALSE, row.names=1)
samples <- c("wt01", "wt02", "wt03")
conditions <- c("wt01", "wt02", "wt03")

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
