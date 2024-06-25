library(tidyverse)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)


# 下载BRCA（乳腺癌）数据

# 查看可下载的癌症类型
# getGDCprojects()$project_id


# 查询数据
query <- GDCquery(
  project = c("TCGA-BRCA"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)


# 下载
GDCdownload(query)

# 读取下载数据
data = GDCprepare(query)

# 保存文件
# saveRDS(data, "data.rds")



# 读取文件
data = readRDS("data.rds")
#提取表达矩阵
data_fpkm <- assay(data, i = 5) %>% 
  as.data.frame()
# 提取注释文件
anno <- rowData(data) 

# 注释
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



# 保存表达矩阵
# write.csv(expr, "expr.csv", row.names = F)


# 提取临床信息文件
pd = colData(data)
# saveRDS(pd,"pd.rds")

