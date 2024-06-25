library(tidyverse)
library(dplyr)
library(ggvenn)


# 
ts = read.csv("targetscan_outfile.csv", header = T)
rdb = read.csv("miRDB target prediction data.csv", header = T)
randa = read.csv("miranda.csv", header = T)

head(ts)
head(rdb)
head(randa)


ts_pair = paste(ts$miRNA_family_ID,"_",ts$a_Gene_ID)
rdb_pair = paste(rdb$miRNA.Name,"_",rdb$Gene.Symbol)
randa_pair = paste(randa$miRNA,"_",randa$gene)

head(ts_pair)
head(rdb_pair)
head(randa_pair)

pair = intersect(intersect(ts_pair,rdb_pair),randa_pair)
pair

data = list("Targetscan" = ts_pair,
                  "miRDB" = rdb_pair,
                  "miRanda" = randa_pair)

ggvenn(data,
       fill_color = c("#3392ff","#f4552f","#55df94"),
       fill_alpha = 0.7)


# write.csv(pair,"pair.csv")


data = read.csv("pair.csv", header = T)
head(data)
dim(data)

head(ts)
ts_select = subset(ts, miRNA_family_ID %in% data$miRNA & a_Gene_ID%in% data$Gene)
dim(ts_select)

head(rdb)
rdb_select = subset(rdb, miRNA.Name %in% data$miRNA & Gene.Symbol%in% data$Gene)
dim(rdb_select)

head(randa)
randa_select = subset(randa, miRNA %in% data$miRNA & gene %in% data$Gene)
dim(randa_select)
