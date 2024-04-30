getwd()

countData <- read.csv("all_featureCountOutput.Rmatrix.txt", header=TRUE, row.names=1, sep="\t")
dim(countData)
head(countData) 

barplot(colSums(countData)*1e-6,
        names=colnames(countData),
        ylab="Library size (millions)")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(pheatmap)

# create experiment labels (two conditions)
colData <- DataFrame(condition=factor(c("YRI", "YRI", "YRI", "GBR", "GBR","YRI", "GBR", "GBR", "GBR", "GBR", "YRI", "YRI")))

dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))

# run DEseq
dds <- DESeq(dds)

plotMA(dds)

# get differentially expressed genes
res <- results(dds)

# order by BH adjusted p-value
resOrdered <- res[order(res$padj),]

head(resOrdered)
nrow(sig)


# how many differentially expressed genes ? FDR=10%, |fold-change|>2 (up and down)
# get differentially expressed gene matrix
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    abs(resOrdered$log2FoldChange)>=1,]

# top of the differentially expressed genes
head(sig)

selected <- rownames(sig);selected

# Checking for any NA, NaN, Inf values
sum(is.na(data_for_heatmap))
sum(is.nan(data_for_heatmap))
sum(is.infinite(data_for_heatmap))

data_for_heatmap <- log2(counts(dds, normalized = TRUE)[rownames(dds) %in% selected,])

# Remove columns with any NA/NaN/Inf
data_for_heatmap <- data_for_heatmap[, colSums(is.na(data_for_heatmap) | is.nan(data_for_heatmap) | is.infinite(data_for_heatmap)) == 0]

# Alternatively, remove rows with any NA/NaN/Inf
data_for_heatmap <- data_for_heatmap[rowSums(is.na(data_for_heatmap) | is.nan(data_for_heatmap) | is.infinite(data_for_heatmap)) == 0, ]


hmcol <- colorRampPalette(c("blue", "white", "red"))(100)

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
