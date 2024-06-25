library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(readxl)
library(cols4all)
library(tidyr)


# 
data = read.csv("TCGA_expr.csv", row.names = 1, header = T)
# 进行log转换
data = log(data+1)


# 
Sample = colnames(data)
length(Sample) # 
tumor = grep("^TCGA.\\w+\\.\\w+.\\.0\\d+\\w", Sample, value = TRUE)
norm = grep("^TCGA.\\w+\\.\\w+.\\.1\\d+\\w", Sample, value = TRUE)
length(norm) # 
length(tumor) # 
# 
group <- c(rep('Normal', length(norm)), rep('Tumor', length(tumor)))
# 
info1 <- data.frame(Sample = c(norm, tumor), Group = group)
head(info1)
row.names(info1) <- info1[,1]
info2 <- info1[,-1]   #

# 
geneset = read.csv("migration_geneset.csv",header = T)
head(geneset)
gene_list = split(geneset$Gene.Symbol,geneset$Phenotype)
head(gene_list)

# ssGSEA
re <- gsva(as.matrix(data), gene_list, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)
# 
# write.csv(re,"ssGSEA_score.csv")

