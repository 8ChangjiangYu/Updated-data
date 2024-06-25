library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)

# 
gene = read.csv("wgcna_gene.csv", header = T)
head(gene)

table(gene$moduleColor)
# blue floralwhite     magenta        pink

gene_blue = subset(gene, moduleColor == "blue")$gene_name
gene_floralwhite = subset(gene, moduleColor == "floralwhite")$gene_name
gene_magenta = subset(gene, moduleColor == "magenta")$gene_name
gene_pink = subset(gene, moduleColor == "pink")$gene_name


# 
bp_go <- function(gene) {
  gene_id = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  enrichGO(gene_id$ENTREZID,
           OrgDb = "org.Hs.eg.db",
           keyType = "ENTREZID",
           ont = "BP",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.2)
}

# 
bp_blue = bp_go(gene_blue)
barplot(bp_blue)

bp_floralwhite = bp_go(gene_floralwhite)
barplot(bp_floralwhite)

bp_magenta = bp_go(gene_magenta)
barplot(bp_magenta)

bp_pink = bp_go(gene_pink)
barplot(bp_pink)


# 
# write.csv(bp_blue, "bp_blue.csv")
# write.csv(bp_floralwhite, "bp_floralwhite.csv")
# write.csv(bp_magenta, "bp_magenta.csv")
# write.csv(bp_pink, "bp_pink.csv")


# 
rm(list = ls())
bp_blue = read.csv("bp_blue_select.csv", header = T, row.names = 1)
bp_magenta = read.csv("bp_magenta_select.csv", header = T, row.names = 1)
bp_pink = read.csv("bp_pink_select.csv", header = T, row.names = 1)

p1 = ggplot(bp_blue, aes(Count,reorder(Description,Count), Count, size = Count, color = pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#f179cd",high = "#97eec6") +
  labs(title = "Biological process of MEblue genes", y = "Description")


p2 = ggplot(bp_magenta, aes(Count,reorder(Description,Count), Count, size = Count, color = pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#f179cd",high = "#97eec6") +
  labs(title = "Biological process of MEmagenta genes", y = "Description")


p3 = ggplot(bp_pink, aes(Count,reorder(Description,Count), Count, size = Count, color = pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#f179cd",high = "#97eec6") +
  labs(title = "Biological process of MEpink genes", y = "Description")
p1;p2;p3

GSEA()
