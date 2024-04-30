setwd("C:/Users/kaaru/Desktop/Assignment/IB2/RNAseq")

data <- read.csv("data.csv", row.names = 1)

install.packages(c("factoextra", "FactoMineR"))

library("factoextra")
library("FactoMineR")

pca.data <- PCA(data[,-1], scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))

fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE) 
pca.data <- PCA(t(data[,-1]), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(pca.data, col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)

install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr) 

a <- fviz_pca_ind(pca.data, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                  repel = TRUE)

ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())

pca.data <- PCA(data[,-1], scale.unit = TRUE,ncp = 2, graph = FALSE)

data$Lineage <- as.factor(data$Lineage)

install.packages("RColorBrewer")
library(RColorBrewer)
nb.cols <- 3
mycolors <- colorRampPalette(brewer.pal(3, "Set1"))(nb.cols)

a <- fviz_pca_ind(pca.data, col.ind = data$Lineage,
                  palette = mycolors, addEllipses = TRUE)

ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cell type", legend.position = "top",
      ggtheme = theme_minimal())
