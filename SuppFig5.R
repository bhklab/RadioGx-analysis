#### Supplementary Figure 5

library(RadioGx)
library(xlsx)
library(gdata)
library(piano)

setwd("PathwayAnalysis-Hypoxia")

load('Oxic_Hypoxia5LQres.RData')
load('CCLERNAseq.RData')

common = intersect(colnames(edata1),rownames(hypoxia.LQ))

ccle.rna1 <- edata1[,common]
Yard1 <- hypoxia.LQ[common,]

# Some samples in the CCLERNAseq data has non zero rows for a given gene
# But when we subset based on the common samples, we are removing the non-zero samples
# and we end up with all zeros for a given gene. So, we need to remove all genes that have zero expression
samp.ind1 <- which(rowSums(ccle.rna1) == 0)
ccle.rna1 <- ccle.rna1[-samp.ind1,]

cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(ccle.rna1)) {
  tmp <- cor.test(ccle.rna1[i,],as.numeric(as.character(Yard1[,'Oxic'])),method="spearman")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(ccle.rna1)
metric <- data.frame(Gene,cor.val,p.val)
metric$Gene <- as.character(metric$Gene)
colnames(metric) <- c("Gene","Correlation","Pvalue")
metric[,4] <- p.adjust(metric$Pvalue,method = "fdr",n = length(metric$Pvalue))
metric.oxic <- subset(metric,metric$V4<0.05)

cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(ccle.rna1)) {
  tmp <- cor.test(ccle.rna1[i,],as.numeric(as.character(Yard1[,'Hyp5'])),method="spearman")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(ccle.rna1)
metric1 <- data.frame(Gene,cor.val,p.val)
metric1$Gene <- as.character(metric1$Gene)
colnames(metric1) <- c("Gene","Correlation","Pvalue")
metric1[,4] <- p.adjust(metric1$Pvalue,method = "fdr",n = length(metric1$Pvalue))
metric1.hyp <- subset(metric1,metric1$V4<0.05)

all.genes <- union(metric.oxic$Gene,metric1.hyp$Gene)
all.genes1 <- data.frame(all.genes, metric$V4[match(all.genes,metric$Gene)],metric1$V4[match(all.genes,metric1$Gene)],NA)
colnames(all.genes1) <- c("Genes","Oxic-FDR","Hyp-FDR","Color")

all.genes1[which(all.genes1[,2] < 0.05 & all.genes1[,3] < 0.05) ,4] <- "red"
all.genes1[which(all.genes1[,2] < 0.05 & all.genes1[,3] > 0.05) ,4] <- "blue"
all.genes1[which(all.genes1[,2] > 0.05 & all.genes1[,3] < 0.05) ,4] <- "magenta"

pdf("OxicHypoxic-FDR.pdf")
plot(-log10(all.genes1[,2]),-log10(all.genes1[,3]),col=all.genes1[,4],xlab="-log10(Oxic-FDR)",ylab="-log10(Hypoxic-FDR)",
     cex.lab=1.5,xlim=c(0.45,9.5),ylim=c(0.45,9.5),cex.axis=1.5)
legend("topleft", c("Significant in Oxic and Hypoxic","Significant in Oxic","Significant in Hypoxic"), 
       xpd = TRUE, bty = "n", pch = c(0,0), col=c("red","blue","magenta"), cex = 1.2)

dev.off()


