library(tidyverse)
library(dplyr)

data = read.delim("mature.fa", header = F)
View(data)

data$V1[seq(2,length(data$V1),by = 2)]

df = data.frame("miRNA" = data$V1[seq(1,length(data$V1),by = 2)],
                "sequence" = data$V1[seq(2,length(data$V1),by = 2)])
df = subset(df, miRNA %in% grep("^>hsa-", df$miRNA, value = T))
head(df)
fasta_content = paste(df$miRNA,"\n",df$sequence)
writeLines(fasta_content, "human_miRNA.fasta")