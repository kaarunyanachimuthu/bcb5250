# working directory
getwd() 

# read in count matrix
countData <- read.csv("count_matrix.txt", header=T, row.names=1, sep="\t") 
dim(countData)
head(countData) 

# basic QC
barplot(colSums(countData)*1e-6,
        names=colnames(countData),
        ylab="Library size (millions)")

# install deseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# load library
library(DESeq2)

# create experiment labels (two conditions)
colData <- DataFrame(condition=factor(c("ctrl","ctrl","ctrl", "treat", "treat", "treat")))

# create DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))

# run DEseq
dds <- DESeq(dds)

# visualize differentially expressed genes
plotMA(dds)

# get differentially expressed genes
res <- results(dds)

# order by BH adjusted p-value
resOrdered <- res[order(res$padj),]

# top of ordered matrix
head(resOrdered)

# how many differentially expressed genes ? FDR=10%, |fold-change|>2 (up and down)
# get differentially expressed gene matrix
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    abs(resOrdered$log2FoldChange)>=1,]

# top of the differentially expressed genes
head(sig)

# how to create a heat map
# select genes
selected <- rownames(sig);selected


data_for_heatmap <- log2(counts(dds, normalized = TRUE)[rownames(dds) %in% selected,])
hmcol <- colorRampPalette(c("blue", "white", "red"))(100)
library(pheatmap)
pheatmap()

pheatmap(data_for_heatmap,
         color = hmcol,
         scale = "row",
         clustering_method = "complete", 
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6, 
         fontsize_col = 12, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         silent_clusters = TRUE 
)
