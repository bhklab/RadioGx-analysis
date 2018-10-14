############################################################
#
# Figure 4    
#
############################################################

# Pathway analysis for alpha/beta ratio in oxic conditions

library(piano)
library(VennDiagram)

setwd("~/RadioGxAnalysis/")

load('ClevelandLQFitVal-GoodCurves.RData')
test <- subset(res2,res2[,3] !=0)
test1 <- test[-which("ZR-75-30" == rownames(test)),] # remove this cell line "ZR-75-30" as the beta value is 10^-19
df <- test1[,2]/test1[,3]
test1 <- cbind(test1,df)
celllines_alphabeta_vals <- test1

# Fig 4A
pdf("Hist_AlphaBetaRatio")
hist(log10(celllines_alphabeta_vals[,6]),col="red",breaks = 50,xlab="log10(alpha/beta)")
dev.off()

# Fig 4B
# Correlation of alpha with AUC, beta with AUC, alpha with beta

e1 <- cor.test(celllines_alphabeta_vals[,2],celllines_alphabeta_vals[,4])
e2 <- cor.test(celllines_alphabeta_vals[,3],celllines_alphabeta_vals[,4])
e3 <- cor.test(celllines_alphabeta_vals[,2],celllines_alphabeta_vals[,3])

# With error bars to the bar plot
F <- as.numeric(c(e1$estimate,e2$estimate,e3$estimate))
L <- as.numeric(c(min(e1$conf.int),min(e2$conf.int),min(e3$conf.int)))
U <- as.numeric(c(max(e1$conf.int),max(e2$conf.int),max(e3$conf.int)))

pdf("CorrPlot-AlphaAUC_BetaAUC_AlphaBeta.pdf",width=8,height=8)
par(oma=c(2,2,2,2)) # bottom,left,top,right
b <- barplot(F, beside=TRUE,names.arg=c("Alpha_AUC","Beta_AUC","Alpha_Beta"),las=0.2,ylab="Correlation (Pearson)",
             col=c("#66C2A5","#FC8D62","#8DA0CB"),cex.axis = 1.75,cex.names = 1.5,cex.lab = 1.75,ylim=c(-1,1))
arrows(x0=b, x1=b, y0=L, y1=U, code=3, angle=90) 
dev.off()

# Fig 4C: Pathway analysis using alpha and beta
load('CCLERNAseq.RData')
common = intersect(colnames(edata1),rownames(celllines_alphabeta_vals))
ccle.rna1 <- edata1[,common]
Yard1 <- celllines_alphabeta_vals[common,]

# Some samples in the CCLERNAseq data has non zero rows for a given gene
# But when we subset based on the common samples, we are removing the non-zero samples
# and we end up with all zeros for a given gene. So, we need to remove all genes that have zero expression
samp.ind1 <- which(rowSums(ccle.rna1) == 0)
ccle.rna1 <- ccle.rna1[-samp.ind1,]

cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(ccle.rna1)) {
  #tmp <- cor.test(ccle.rna1[i,],Yard1[,'LQ-Alpha'],method="spearman")
  tmp <- cor.test(ccle.rna1[i,],Yard1[,'LQ-Beta'],method="spearman")
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
res.gsea <- subset(gseares1,gseares1[,"FDR"] < 0.15) 
save(gsares,file="GSEAResult-Alpha.RData")
#save(gsares,file="GSEAResult-Beta.RData")

################################################################
#
# Find common pathways between AUC, alpha, beta, alpha/beta
#
################################################################

load('GSEAResult-Alpha.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
#write.csv(gseares1,file="Pathways-Alpha-NoCutOff.csv")
res.alpha <- subset(gseares1,gseares1[,"FDR"] < 0.05) 

load('GSEAResult-Beta.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
#write.csv(gseares1,file="Pathways-Beta-NoCutOff.csv")
res.beta <- subset(gseares1,gseares1[,"FDR"] < 0.05) 

load("GSEA-AUC.RData")
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1AUC <- gseares1[,-(4:5)]
res.auc <- subset(gseares1AUC,gseares1AUC[,"FDR"] < 0.05) 

# Venn diagram (alpha, beta, AUC)
pdf("Venn-Alpha_Beta_AUC.pdf")
venn.plot <- draw.triple.venn(
  area1 = nrow(res.alpha),
  area2 = nrow(res.beta),
  area3 = nrow(res.auc),
  n12 = length(intersect(res.alpha$Name,res.beta$Name)),
  n23 = length(intersect(res.beta$Name,res.auc$Name)),
  n13 = length(intersect(res.alpha$Name,res.auc$Name)),
  n123 = length(Reduce(intersect, list(res.alpha$Name,res.beta$Name,res.auc$Name))),
  category = c("Alpha", "Beta", "AUC"),
  fill = c("orange", "red", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red","blue")
)
dev.off()
grid.draw(venn.plot)


