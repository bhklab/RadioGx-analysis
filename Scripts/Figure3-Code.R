##########################################################
# Compare pathways between SF2 and AUC
# Expression data: CCLE RNA seq
# Response variables: SF2 and AUC
# Method: GSEA using piano
#
##########################################################

library(xlsx)
library(gdata)
library(genefu)
library(piano)
library(RColorBrewer)
library(tidyverse)
library(plyr)
require(VennDiagram)

setwd("~/RadioGxAnalysis/")

# Figure 3A

load('DoseResponse.RData')
load('ClevelandLQFitVal-GoodCurves.RData')
Yard <- dose.response[rownames(res2),]

a <- cor.test(Yard$`2Gy`,Yard$YardAUC)

# color the scatter plot of AUC and SF2 by tertiles
b <- data.frame(Yard$`2Gy`,Yard$YardAUC)
colnames(b) <- c("SF2","AUC")
rownames(b) <- rownames(Yard)
b[,3] <- NA
colnames(b)[3] <- "Color"

# Tertiles by range of SF2
#test <- cut(b$SF2, breaks = 3)
#table(test)

b$Color[which(b$SF2>=0.064 & b$SF2<=0.45)] <- "#66C2A5"
b$Color[which(b$SF2>0.45 & b$SF2<=0.836)] <- "#FC8D62"
b$Color[which(b$SF2>0.836 & b$SF2<=1.23)] <- "#8DA0CB"

pdf("SF2-AUC-Tertile.pdf")
cols <- b$Color
plot(b$SF2,b$AUC,cex.lab = 2,cex.axis=2,xlab="SF2",ylab="AUC",pch=16,col = cols)
legend("topleft", c("Sensitive","Intermediate","Resistant"), pch=15, col=c("#66C2A5","#FC8D62","#8DA0CB"), bty="n",cex=1.5)
legend("bottomright", legend=sprintf("Correlation= %.2f (p<e-16)",a$estimate),bty="n",cex=1.5,horiz = TRUE)
dev.off()

############################################################
# Figure 3B

load('DoseResponse.RData')
load('ClevelandLQFitVal-GoodCurves.RData')
Yard <- dose.response[rownames(res2),]

b <- data.frame(Yard$`2Gy`,Yard$YardAUC)
colnames(b) <- c("SF2","AUC")
rownames(b) <- rownames(Yard)

c <- b %>%mutate(quantile = ntile(b$SF2, 3))

# Tertiles by range
#test <- cut(b$SF2, breaks = 3)
#table(test)

c1 <- subset(c,c$SF2>=0.064 & c$SF2<=0.45)
c2 <- subset(c,c$SF2>0.45 & c$SF2<=0.836)
c3 <- subset(c,c$SF2>0.836 & c$SF2<=1.23)

e1 <- cor.test(c1$SF2,c1$AUC)
e2 <- cor.test(c2$SF2,c2$AUC)
e3 <- cor.test(c3$SF2,c3$AUC)

# With error bars to the bar plot
F <- as.numeric(c(e1$estimate,e2$estimate,e3$estimate))
L <- as.numeric(c(min(e1$conf.int),min(e2$conf.int),min(e3$conf.int)))
U <- as.numeric(c(max(e1$conf.int),max(e2$conf.int),max(e3$conf.int)))

pdf("CorrPlot-SF2-AUC-Tertiles.pdf",width=8,height=8)
par(oma=c(2,2,2,2)) # bottom,left,top,right
b <- barplot(F, beside=TRUE,names.arg=c("Sensitive","Intermediate","Resistant"),las=0.2,ylab="Correlation (Pearson)",
             col=c("#66C2A5","#FC8D62","#8DA0CB"),ylim=c(0,1),cex.axis = 1.75,cex.names = 1.5,cex.lab = 1.75)
arrows(x0=b, x1=b, y0=L, y1=U, code=3, angle=90) 
dev.off()

###############################################################
###############################################################

load('DoseResponse.RData')
load('ClevelandLQFitVal-GoodCurves.RData')
Yard <- dose.response[rownames(res2),]

load('CCLERNAseq.RData')
ccle.rna <- edata1

common <- intersect(rownames(Yard),colnames(ccle.rna))
ccle.rna1 <- ccle.rna[,common]
Yard1 <- Yard[common,]

# Some samples in the CCLERNAseq data has non zero rows for a given gene
# But when we subset based on the common samples, we are removing the non-zero samples
# and we end up with all zeros for a given gene. So, we need to remove all genes that have zero expression
samp.ind1 <- which(rowSums(ccle.rna1) == 0)
ccle.rna1 <- ccle.rna1[-samp.ind1,]

# Pathway analysis using SF2 values
cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(ccle.rna1)) {
  tmp <- cor.test(ccle.rna1[i,],as.numeric(as.character(Yard1$`2Gy`)),method="spearman")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(ccle.rna1)
metric <- data.frame(Gene,cor.val,p.val)
metric$Gene <- as.character(metric$Gene)
colnames(metric) <- c("Gene","Correlation","Pvalue")

load('GeneSets-IPA-FinalVersion-EntID.RData')
a <- loadGSC(gSets_IPA_EntID)

stats <- as.vector(metric$Correlation)
names(stats) <- as.vector(metric$Gene)
gsares <- runGSA(geneLevelStats=stats,gsc=a,nPerm=10000,geneSetStat="gsea",adjMethod="none")
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
res.gsea <- subset(gseares1,gseares1[,"FDR"] < 0.05) 
save(gsares,file="GSEA-AUC.RData")

# Pathway analysis using Yard AUC values
cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(ccle.rna1)) {
  tmp <- cor.test(ccle.rna1[i,],as.numeric(as.character(Yard1$YardAUC)),method="spearman")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(ccle.rna1)
metric <- data.frame(Gene,cor.val,p.val)
metric$Gene <- as.character(metric$Gene)
colnames(metric) <- c("Gene","Correlation","Pvalue")

load('GeneSets-IPA-FinalVersion-EntID.RData')
a <- loadGSC(gSets_IPA_EntID)

stats <- as.vector(metric$Correlation)
names(stats) <- as.vector(metric$Gene)
gsares <- runGSA(geneLevelStats=stats,gsc=a,nPerm=10000,geneSetStat="gsea",adjMethod="none")
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
res.gsea <- subset(gseares1,gseares1[,"FDR"] < 0.05) 
save(gsares,file="GSEA-SF2.RData")
################################################################

# Fig3C:venn diagram

pdf("Venn_SF2-AUC.pdf")
draw.pairwise.venn(area1 = nrow(res.auc), area2 = nrow(res.sf2), cross.area = length(intersect(res.auc$Name,res.sf2$Name)),
                   category = c("AUC","SF2"),lty = rep("blank", 2), fill = c("red","blue"))
dev.off()


#######################################################
# Figure 3D: Plot -log10FDRAUC vs. -log1-FDRSF2 by SF2
# color by common and diff pathways

load('GSEA-AUC.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1AUC <- gseares1[,-(4:5)]
res.auc <- subset(gseares1AUC,gseares1AUC[,"FDR"] < 0.05) 

load('GSEA-SF2.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1SF2 <- gseares1[,-(4:5)]
res.sf2 <- subset(gseares1SF2,gseares1SF2[,"FDR"] < 0.05) 

all.pathways <- matrix(NA,length(union(res.auc$Name,res.sf2$Name)),3)
rownames(all.pathways) <- union(res.auc$Name,res.sf2$Name)
colnames(all.pathways) <- c("AUC-FDR","SF2FDR","Color")
all.pathways[,1] <- gseares1AUC$FDR[match(rownames(all.pathways),gseares1AUC$Name)]
all.pathways[,2] <- gseares1SF2$FDR[match(rownames(all.pathways),gseares1SF2$Name)]

common <- intersect(res.auc$Name,res.sf2$Name)
all.pathways[match(common,rownames(all.pathways)),3] <- "#8DA0CB"

diff.path.auc <- setdiff(res.auc$Name,res.sf2$Name)
diff.path.sf2 <- setdiff(res.sf2$Name,res.auc$Name)

all.pathways[match(diff.path.auc,rownames(all.pathways)),3] <- "red"
all.pathways[match(diff.path.sf2,rownames(all.pathways)),3] <- "blue"

a <- data.frame(-log10(as.numeric(all.pathways[,1])),-log10(as.numeric(all.pathways[,2])),all.pathways[,3])
rownames(a) <- rownames(all.pathways)
colnames(a) <- c("AUC","SF2","Color")
a$Color <- as.character(a$Color)

pdf("ComparePathways_SF2AUC.pdf")
plot(a$AUC,a$SF2,col=a$Color,pch=19,xlab="-log10(FDR_AUC)",ylab="-log10(FDR_SF2)",xlim=c(0.4,2.2),ylim=c(0.4,2.2),
     cex.lab=1.5,cex.axis = 1.5)
legend("topleft",fill = c("#8DA0CB","red","blue"),cex=1.2,
       legend = c("Common pathways (35)","Pathways specific to AUC (42)","Pathways specific to SF2 (3)"),bty = "n")
dev.off()










