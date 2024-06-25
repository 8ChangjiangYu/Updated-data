# 
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)

# 
expr = read.csv("tcga_expr.csv", header = T, row.names = 1)
expr = log2(expr+1)

# 
pd = read.delim("TCGA-BRCA.survival.tsv")
head(pd)
pd = select(pd, c("sample","OS","OS.time"))
colnames(pd) = c("sample","status","time")
head(pd)

# 
hub_gene = read.csv("hubgene.csv", header = T)
head(hub_gene)


# 
expr = expr %>% 
  t() %>% 
  as.data.frame()
hub_expr = select(expr, hub_gene$Name)
head(hub_expr)


# 
dim(pd) # 1260    3
head(pd)
dim(hub_expr) # 1231   40
hub_expr$sample = substr(rownames(hub_expr),1,16)
pd$sample = gsub("-",".",pd$sample)


data = merge(hub_expr, pd, by = "sample")
dim(data) # 1208   43
head(data)

# write.csv(data, "surrvive_data.csv", row.names = F)




# 
rm(list = ls())
data = read.csv("surrvive_data.csv", header = T)
head(data)

dim(data) # 1208   43
for (i in colnames(data)[2:41]){
  name = paste0(i,"_","group")
  # value = data[,i]
  data[, name] = ifelse(data[,i] > median(data[,i]), "high", "low")
}
head(data)

# 
data_value = data.frame()
for (i in colnames(data)[2:41]){
  name = paste0(i,"_","group")
  formula_str <- paste("Surv(time, status) ~", name)
  cox_model1 <- coxph(as.formula(formula_str), data = data)
  data_value[i,"pvalue"] = summary(cox_model1)$coefficients[,5]
}
head(data_value)
dim(data_value) # 
data_value_diff = filter(data_value, pvalue < 0.05)
dim(data_value_diff) # 

# 
gene = rownames(data_value_diff)
gene # "MRPL20"   "MRPL12"   "AURKAIP1" "NDUFB7"   "ATP5F1D"  "COL4A1"   "BGN"      "VWF" 

head(data)
## 
# 
df = data %>% select(c(gene, "MRPL20_group","MRPL12_group","AURKAIP1_group","NDUFB7_group",
                    "ATP5F1D_group","COL4A1_group","BGN_group","VWF_group","status","time"))
head(df)
dim(df)

for (i in colnames(df)[1:8]){
  name = paste0(i,"_","group")
  formula_str <- paste("Surv(time, status) ~", name)
  cox_model2 <- survfit(as.formula(formula_str), data = df)
  title = paste0(i," Survival Curve")
  p = ggsurvplot(cox_model2, # 
             data = df,  # 
             conf.int = TRUE, # 
             pval = TRUE, # 
             risk.table = TRUE, # 
             surv.median.line = "hv", # 
             # add.all = TRUE, # 
             palette = "hue") +  #  
    ggtitle(title)
  print(p)
}


cox_model1 = survfit(Surv(time, status) ~ MRPL20_group, data = df)
cox_model2 = survfit(Surv(time, status) ~ MRPL12_group, data = df)
cox_model3 = survfit(Surv(time, status) ~ AURKAIP1_group, data = df)
cox_model4 = survfit(Surv(time, status) ~ NDUFB7_group, data = df)
cox_model5 = survfit(Surv(time, status) ~ ATP5F1D_group, data = df)
cox_model6 = survfit(Surv(time, status) ~ COL4A1_group, data = df)
cox_model7 = survfit(Surv(time, status) ~ BGN_group, data = df)
cox_model8 = survfit(Surv(time, status) ~ VWF_group, data = df)



p1 = ggsurvplot(cox_model1, # 
           data = df,  # 
           conf.int = TRUE, # 
           pval = TRUE, # 
           risk.table = TRUE, # 
           surv.median.line = "hv", # 
           # add.all = TRUE, # 
           palette = "hue") +  #  
  ggtitle("MRPL20 Survival Curve")

p2 = ggsurvplot(cox_model2, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  # 
  ggtitle("MRPL12 Survival Curve")

p3 = ggsurvplot(cox_model3, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  #  
  ggtitle("AURKAIP1 Survival Curve")

p4 = ggsurvplot(cox_model4, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  #  
  ggtitle("NDUFB7 Survival Curve")

p5 = ggsurvplot(cox_model5, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  # 
  ggtitle("ATP5F1D Survival Curve")

p6 = ggsurvplot(cox_model6, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  #  
  ggtitle("COL4A1 Survival Curve")

p7 = ggsurvplot(cox_model7, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  # 
  ggtitle("BGN Survival Curve")

p8 = ggsurvplot(cox_model8, # 
                data = df,  # 
                conf.int = TRUE, # 
                pval = TRUE, # 
                risk.table = TRUE, # 
                surv.median.line = "hv", # 
                # add.all = TRUE, # 
                palette = "hue") +  # 
  ggtitle("VWF Survival Curve")

