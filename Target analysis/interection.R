library(tidyverse)
library(dplyr)

data = read.csv("DGIdb_gene_interaction_results-2023_12_29 09_24_48.tsv",sep = "\t")
View(data)


ggplot(data, aes(gene,drug,size = interaction.score, color = regulatory.approval))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("#D1392B","#A9C4E6"))
# write.csv(data,"Table1.csv", row.names = F)



