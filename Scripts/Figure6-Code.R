##########################################################
#
# Figure 6: Tissue specificity of radiation response
# GSEA using IPA defined pathways using PIANO package
#
##########################################################

library(xlsx)
library(genefu)
library(piano)
library(gplots)
library(RColorBrewer)

setwd("~/RadioGxAnalysis/")

load('ClevelandLQFitVal-GoodCurves.RData')
load('YardHistology.RData')

Histology <- YardHistology[match(rownames(res2),rownames(YardHistology)),]
Histology$Primarysite[which(Histology$Primarysite == "")] <- "Others"

tissue <- data.frame(table(Histology$Primarysite))
tissue <- subset(tissue,tissue$Freq>15)
tissue$Var1 <- as.character(tissue$Var1)

load('CCLERNAseq.RData')

for (k in 1:nrow(tissue)) {
  print(k)
  cell.lines <- rownames(Histology)[which(tissue$Var1[k] == Histology$Primarysite)]
  #tmp <- match(cell.lines,colnames(edata1)) # some cell lines in Yard do not have gene expression in the CCLE
  cell.lines <- cell.lines[which(is.na(match(cell.lines,colnames(edata1))) == F)] # some cell lines in Yard do not have gene expression in the CCLE
  tissue.edata <- edata1[,cell.lines]
  
  # Some samples in the CCLERNAseq data has non zero rows for a given gene
  # But when we subset based on the common samples, we are removing the non-zero samples
  # and we end up with all zeros for a given gene. So, we need to remove all genes that have zero expression
  samp.ind1 <- which(rowSums(tissue.edata) == 0)
  tissue.edata1 <- tissue.edata[-samp.ind1,]
  
  auc.values <- res2[cell.lines,'LQ-AUC']
  
  cor.val <- as.numeric()
  p.val <- as.numeric()
  for (i in 1:nrow(tissue.edata1)) {
    tmp <- cor.test(tissue.edata1[i,],as.numeric(as.character(auc.values)),method="spearman")
    cor.val[i] <- as.numeric(tmp$estimate)
    p.val[i] <- as.numeric(tmp$p.value)
  }
  Gene <- rownames(tissue.edata1)
  metric <- data.frame(Gene,cor.val,p.val)
  metric$Gene <- as.character(metric$Gene)
  colnames(metric) <- c("Gene","Correlation","Pvalue")
  
  load('GeneSets-IPA-FinalVersion-EntID.RData')
  a <- loadGSC(gSets_IPA_EntID)
  
  stats <- as.vector(metric$Correlation)
  names(stats) <- as.vector(metric$Gene)
  gsares <- runGSA(geneLevelStats=stats,gsc=a,nPerm=10000,geneSetStat="gsea",adjMethod="none")
  save(gsares,file=paste("GSEA",tissue$Var1[k],".RData",sep=""))
}

#### GSEA on all tissues 

load('YardHistology.RData')
Histology <- YardHistology[match(rownames(res2),rownames(YardHistology)),]
Histology$Primarysite[which(Histology$Primarysite == "")] <- "Others"

others.celllines <- rownames(Histology)[which(Histology$Primarysite == "Others")]
Histology <- Histology[-match(others.celllines,rownames(Histology)),]

load('ClevelandLQFitVal-GoodCurves.RData')
res2 <- res2[-match(others.celllines,rownames(res2)),]  

load('CCLERNAseq.RData')
common.cell.lines <- intersect(rownames(res2),colnames(edata1)) # some cell lines in Yard do not have gene expression in the CCLE

res2 <- res2[match(common.cell.lines,rownames(res2)),] # some cell lines in Yard do not have gene expression in the CCLE
tissue.edata <- edata1[,common.cell.lines]

# Some samples in the CCLERNAseq data has non zero rows for a given gene
# But when we subset based on the common samples, we are removing the non-zero samples
# and we end up with all zeros for a given gene. So, we need to remove all genes that have zero expression
samp.ind1 <- which(rowSums(tissue.edata) == 0)
tissue.edata1 <- tissue.edata[-samp.ind1,]
auc.values <- res2[common.cell.lines,'LQ-AUC']

cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(tissue.edata1)) {
  tmp <- cor.test(tissue.edata1[i,],as.numeric(as.character(auc.values)),method="spearman")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(tissue.edata1)
metric <- data.frame(Gene,cor.val,p.val)
metric$Gene <- as.character(metric$Gene)
colnames(metric) <- c("Gene","Correlation","Pvalue")

load('GeneSets-IPA-FinalVersion-EntID.RData')
a <- loadGSC(gSets_IPA_EntID)

stats <- as.vector(metric$Correlation)
names(stats) <- as.vector(metric$Gene)
gsares <- runGSA(geneLevelStats=stats,gsc=a,nPerm=10000,geneSetStat="gsea",adjMethod="none")
save(gsares,file="GSEA-AllTissues.RData")

#####################################################
# Fig 6A

################################################################

# GSEA of 12 tissues with IPA pathways
# Pathways of 12 tissues 
# FDR < 0.05

################################################################

#breast
load('GSEAbreast.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
breast.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#cns
load('GSEAcentral_nervous_system.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
cns.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#endo
load('GSEAendometrium.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
endo.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#largeint
load('GSEAlarge_intestine.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
largeint.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#oeso
load('GSEAoesophagus.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
oeso.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#liver
load('GSEAliver.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
liver.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#ovary
load('GSEAovary.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
ovary.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#lung
load('GSEAlung.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
lung.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#skin
load('GSEAskin.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
skin.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#pancreas
load('GSEApancreas.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
panc.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#UpAeroTract
load('GSEAupper_aerodigestive_tract.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
uptract.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#UrTract
load('GSEAurinary_tract.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
urtract.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#Stomach
load('GSEAstomach.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
stomach.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#All Tissues
load('GSEA-AllTissues.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
AllTissues.pathways <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 184 pathways are significantly enriched

#####################
all.pathways <- Reduce(union, list(breast.pathways$Name, cns.pathways$Name, largeint.pathways$Name,oeso.pathways$Name,endo.pathways$Name,
                                   liver.pathways$Name,ovary.pathways$Name,lung.pathways$Name,skin.pathways$Name,
                                   panc.pathways$Name,uptract.pathways$Name,urtract.pathways$Name,stomach.pathways$Name,AllTissues.pathways$Name))
mat <- matrix(NA,length(all.pathways),14)
colnames(mat) <- c("Breast","CNS","LargeInt","Endometrium","Oesophagus","Liver","Ovary","Lung","Skin",
                   "Pancreas","UpAero","UrTract","Stomach","AllTissues")
rownames(mat) <- all.pathways

mat[,"Breast"] <- breast.pathways$'Stat (dist.dir)'[match(rownames(mat),breast.pathways$Name)]
mat[,"CNS"] <- cns.pathways$'Stat (dist.dir)'[match(rownames(mat),cns.pathways$Name)]
mat[,"LargeInt"] <- largeint.pathways$'Stat (dist.dir)'[match(rownames(mat),largeint.pathways$Name)]
mat[,"Endometrium"] <- endo.pathways$'Stat (dist.dir)'[match(rownames(mat),endo.pathways$Name)]
mat[,"Oesophagus"] <- oeso.pathways$'Stat (dist.dir)'[match(rownames(mat),oeso.pathways$Name)]
mat[,"Liver"] <- liver.pathways$'Stat (dist.dir)'[match(rownames(mat),liver.pathways$Name)]
mat[,"Ovary"] <- ovary.pathways$'Stat (dist.dir)'[match(rownames(mat),ovary.pathways$Name)]
mat[,"Lung"] <- lung.pathways$'Stat (dist.dir)'[match(rownames(mat),lung.pathways$Name)]
mat[,"Skin"] <- skin.pathways$'Stat (dist.dir)'[match(rownames(mat),skin.pathways$Name)]
mat[,"Pancreas"] <- panc.pathways$'Stat (dist.dir)'[match(rownames(mat),panc.pathways$Name)]
mat[,"UpAero"] <- uptract.pathways$'Stat (dist.dir)'[match(rownames(mat),uptract.pathways$Name)]
mat[,"UrTract"] <- urtract.pathways$'Stat (dist.dir)'[match(rownames(mat),urtract.pathways$Name)]
mat[,"Stomach"] <- stomach.pathways$'Stat (dist.dir)'[match(rownames(mat),stomach.pathways$Name)]
mat[,"AllTissues"] <- AllTissues.pathways$'Stat (dist.dir)'[match(rownames(mat),AllTissues.pathways$Name)]

colwithNA <- which(apply(mat, 2, function(x) all(is.na(x))))
mat1 <- mat[,-colwithNA]

# Heatmap
mat <- read.csv('Tissue Specific with categories.csv')
mat <- mat[,-16]

mat.all <- read.csv('AllTissues.Pathways.csv')
mat.all <- mat.all[,c(1,2,4)]
colnames(mat.all)[3] <- "AllTissues"

all.pathways <- Reduce(union, list(mat$Name,mat.all$Name))

categories <- matrix(NA,length(all.pathways),2)
colnames(categories) <- c("Category","Name")
categories[,"Name"] <- all.pathways
categories[,"Category"] <- as.character(mat$Category[as.numeric(match(categories[,"Name"],mat$Name))])
categories[282:294,"Category"] <- as.character(mat.all$Category[as.numeric(match(categories[282:294,"Name"],mat.all$Name))])

mat_test <- matrix(NA,length(all.pathways),16)
colnames(mat_test) <- c("Category","Name","Breast","CNS","LargeInt","Endometrium","Oesophagus","Liver","Ovary","Lung","Skin",
                        "Pancreas","UpAero","UrTract","Stomach","AllTissues")
rownames(mat_test) <- all.pathways

mat_test[,"Category"] <- categories[,1]
mat_test[,"Name"] <- all.pathways
mat_test[,"Breast"] <- mat$Breast[match(rownames(mat_test),mat$Name)]
mat_test[,"CNS"] <- mat$CNS[match(rownames(mat_test),mat$Name)]
mat_test[,"LargeInt"] <- mat$LargeInt[match(rownames(mat_test),mat$Name)]
mat_test[,"Endometrium"] <- mat$Endometrium[match(rownames(mat_test),mat$Name)]
mat_test[,"Oesophagus"] <- mat$Oesophagus[match(rownames(mat_test),mat$Name)]
mat_test[,"Liver"] <- mat$Liver[match(rownames(mat_test),mat$Name)]
mat_test[,"Ovary"] <- mat$Ovary[match(rownames(mat_test),mat$Name)]
mat_test[,"Lung"] <- mat$Lung[match(rownames(mat_test),mat$Name)]
mat_test[,"Skin"] <- mat$Skin[match(rownames(mat_test),mat$Name)]
mat_test[,"Pancreas"] <- mat$Pancreas[match(rownames(mat_test),mat$Name)]
mat_test[,"UpAero"] <- mat$UpAero[match(rownames(mat_test),mat$Name)]
mat_test[,"UrTract"] <- mat$UrTract[match(rownames(mat_test),mat$Name)]
mat_test[,"Stomach"] <- mat$Stomach[match(rownames(mat_test),mat$Name)]
mat_test[,"AllTissues"] <- mat.all$AllTissues[match(rownames(mat_test),mat.all$Name)]

biol <- levels(factor(categories[,"Category"]))

tmp <- list()
for (i in 1:length(biol)) {
  tmp[[i]] <- mat_test[which(biol[i] == mat_test[,"Category"]),]
}
res <- do.call("rbind", tmp)
rownames(res) <- as.character(res[,"Name"])
res1 <- res[,-c(1,2,15)]

col = brewer.pal(12, "Set3") 
#col1 <- colorRampPalette(col)(25)
n <- 17 #number of categories
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col1 <- col_vector

mymat <- res1
mymat <- as.matrix(mymat)
mymat[mymat == 0] <- NA
mymat1 <- apply(mymat, 2, as.numeric)

mydf <- data.frame(as.character(res[,"Name"]),as.character(res[,"Category"]))
colnames(mydf) <- c("gene","category")

breaks=0:294
mycol <- colorpanel(n=length(breaks)-1,low="darkblue",mid="grey50",high="yellow")
pdf("RadRes2.pdf", height=10, width=15)
par(mar=c(10, 10, 4, 2) + 0.1)
heatmap.2(mymat1, col=mycol,
          trace="none",
          keysize=1,labRow = FALSE,
          #margins=c(12,5),
          scale="none",
          dendrogram="column",
          Rowv = FALSE,
          #main="Pathways grouped by categories",
          RowSideColors=col1[as.numeric(mydf$category)],na.color = "gray90",
          density.info = "none",cexRow=5,sepwidth=c(0.0000001,0.0000001),
          sepcolor="black",colsep=1:ncol(mymat1))
#legend("bottomleft",legend = biol,col = col1,lty= 1,lwd = 8,cex=0.5)
legend(y=0.75, x=-.15, xpd=TRUE,legend = biol,col = col1, lty= 1,lwd = 5,cex=.7)
dev.off()

######################################################################
# Fig 6B
# Tissue specificity of log10(alphabeta)

load('ClevelandLQFitVal-GoodCurves.RData')
test <- subset(res2,res2[,3] !=0)
test1 <- test[-which("ZR-75-30" == rownames(test)),] # remove this cell line "ZR-75-30" as the beta value is 10^-19
df <- test1[,2]/test1[,3]
test1 <- cbind(test1,df)

celllines_alphabeta_vals <- test1
colnames(celllines_alphabeta_vals)[6] <- c("AlphaBetaRatio")

load('YardHistology.RData')
histlo <- YardHistology$Primarysite[match(rownames(celllines_alphabeta_vals),rownames(YardHistology))]

celllines_alphabeta_vals1 <- data.frame(celllines_alphabeta_vals[,6],histlo)
colnames(celllines_alphabeta_vals1) <- c("AphaBetaRatio","histlo")
celllines_alphabeta_vals1$histlo <- as.character(celllines_alphabeta_vals1$histlo) 

df <- celllines_alphabeta_vals1
tt <- table(df$histlo)
df2 <- subset(df, df$histlo %in% names(tt[tt > 5]))
df3 <- df2[-which(df2$histlo == ""),]

df4 <- subset(df3,df3$AphaBetaRatio != 0)
df4$AphaBetaRatio <- log10(df4$AphaBetaRatio)

long <- melt(df4, id.vars = "histlo")
bymedian <- with(df4, reorder(df4$histlo, df4$AphaBetaRatio, median)) # order by median of oxic
long$histlo <- factor(bymedian)

labels1 <- c("Ovary","Oesophagus","UpAero","Skin","Pancreas","CNS","Liver","Lung",
             "UrTract","Thyroid","Breast","LargeInt","Endometrium")

pdf("AlphaBetaRatioHistology.pdf",width=10,height=10)
par(oma=c(15,10,5,5)) # bottom,left,top,right
boxplot(value ~ variable + histlo, data = long,las=2,col=rep(c("#66C2A5"),length(unique(df4$histlo))),
        at = seq(1,13),xaxt='n',ylab="log10(Alpha/Beta)",
        cex.names = 1.75,cex.axis = 1.75,cex.lab=1.75)
axis(1,at = seq(1,13),labels = labels1,las=2,cex.axis = 1.75)
dev.off()












