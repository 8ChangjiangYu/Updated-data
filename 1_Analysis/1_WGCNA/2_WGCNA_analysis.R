rm(list = ls())

# 
library(dplyr)
library(WGCNA)
options(stringsAsFactors = F)

# 
expr = read.csv("tcga_expr.csv", header = T, row.names = 1)
expr = log2(expr+1)
expr = expr %>% 
  t() %>% 
  as.data.frame()

# 
Sample = rownames(expr)
length(Sample) # 
tumor = grep("^TCGA.\\w+\\.\\w+.\\.0\\d+\\w", Sample, value = TRUE)
length(tumor) # 
expr1 = subset(expr, rownames(expr) %in% tumor)
dim(expr1) # 1118 19938

# 
gsg = goodSamplesGenes(expr1, verbose = 3)
gsg$allOK # http://127.0.0.1:34275/graphics/plot_zoom_png?width=1920&height=1017

# 
if(!gsg$allOK){
  if(sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(expr1)[!gsg$goodGenes], collapse = ", ")));
  if(sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expr1)[!gsg$goodSamples], collapse = ", ")))
  expr = expr[gsg$goodSamples, gsg$goodGenes]
}

# 
sampleTree = hclust(dist(expr1), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)




# 
ssgsea = read.csv("ssGSEA_score.csv", header = T, row.names = 1)
ssgsea = ssgsea %>% 
  t() %>% 
  as.data.frame()
head(ssgsea)
ssgsea1 = subset(ssgsea, rownames(ssgsea) %in% tumor)
dim(ssgsea1) # 1118    4


# 
sampleTree2 = hclust(dist(expr1), method = "average")
traitColors = numbers2colors(ssgsea, signed = FALSE) # 
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels = names(expr1),
                    main = "Sample dendrogram and trait heatmap")
save(expr1, ssgsea, file = "wgcna_datainput.Rdata")


# 
load(file = "wgcna_datainput.Rdata")

# 
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(expr1,powerVector = powers, verbose = 5)
sft

# 
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9 #字符大小

# 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red") #

# 

# 
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")

# 
sft$powerEstimate # 5


# 
cor = WGCNA::cor #
net = blockwiseModules(
  expr1, power = 5,
  TOMType = "unsigned",minModuleSize = 30,
  reassignThreshold = 0,mergeCutHeight = 0.25,
  numericLabels = TRUE,pamRespectsDendro = FALSE,
  saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM",
  maxBlockSize = 5000,
  verbose = 3) # 
cor = stats::cor #



# 
table(net$colors)

# 
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],net$dendrograms[[1]],
                    "Module colors",
                    dendroLabels = FALSE,hang = 0.03,
                    addGuide = TRUE,guideHang = 0.05)


# 
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "modulate.Rdata")
dim(MEs)
# 



# 
# 
load(file = "wgcna_datainput.Rdata")
load(file = "modulate.Rdata")


# 
nGenes = ncol(expr1);
nSamples = nrow(expr1);
# 
MEs0 = moduleEigengenes(expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, ssgsea1, use = "p")#
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 
sizeGrWindow(10,6)
# 
textMatrix = paste(signif(moduleTraitCor, 2),"\n(",
                   signif(moduleTraitPvalue, 1),")",sep = "");
dim(textMatrix)=dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(ssgsea1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))




# 
head(ssgsea1)
cell_migration = as.data.frame(ssgsea1$`increased directional cell migration`)
names(cell_migration) = "cell migration"
head(cell_migration)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(expr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(expr1,cell_migration,use = "p"))#
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(cell_migration), sep = "");
names(GSPvalue) = paste("p.GS.", names(cell_migration),sep = "");
# MMPvalue
View(MMPvalue)
# GSPvalue
View(GSPvalue)



# 
# 
module = c("magenta", "pink", "floralwhite", "blue")
column = match(module, modNames);
moduleGenes = moduleColors %in% module

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in",module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)

names(expr1)#
names(expr1)[moduleColors %in% module] #
geneInfo = data.frame(gene_name = names(expr1)[moduleColors %in% module],
                      moduleColor = moduleColors[moduleColors %in% module])
head(geneInfo)
table(geneInfo$moduleColor)

# 
# write.table(geneInfo, "wgcna_gene.txt", sep = "\t", row.names = F, quote = F)
