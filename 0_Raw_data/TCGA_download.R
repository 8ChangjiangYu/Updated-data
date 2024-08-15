library(tidyverse)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)


# 

# 
# getGDCprojects()$project_id


# 
query <- GDCquery(
  project = c("TCGA-BRCA"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)


# 
GDCdownload(query)

# 
data = GDCprepare(query)

# 
# saveRDS(data, "data.rds")



# 
data = readRDS("data.rds")
#
data_fpkm <- assay(data, i = 5) %>% 
  as.data.frame()
# 
anno <- rowData(data) 

# 
anno_info = anno %>% 
  as.data.frame() %>% 
  subset(gene_type == "protein_coding") %>% 
  select(gene_id, gene_name)
head(anno_info)

data_fpkm$gene_id = rownames(data_fpkm)
expr = merge(data_fpkm, anno_info, by = "gene_id")
expr = expr %>% 
  select(-gene_id) %>% 
  select(gene_name, everything())
expr[1:3,1:5]

length(unique(str_sub(colnames(expr),1,16)))



# 
# write.csv(expr, "expr.csv", row.names = F)


# 
pd = colData(data)
# saveRDS(pd,"pd.rds")

