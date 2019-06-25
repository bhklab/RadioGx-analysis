##########################################################
#
# Figure 8: RT and Chemo genetic dependencies using PCL
#
##########################################################

library(PharmacoGx)
library(RadioGx)
library(xlsx)
library(gdata)

load('CTRPv2AUCPublishedVals.RData')
load('DoseResponse.RData')
load('ClevelandLQFitVal-GoodCurves.RData')

dose.response <- dose.response[rownames(res2),]

common <- intersect(rownames(dose.response),colnames(CTRPv2auc))
Yard <- dose.response[common,]
CTRPv2auc1 <- t(CTRPv2auc[,common])

z1.cor <- NULL
z1.pval <- NULL

for(i in 1:ncol(CTRPv2auc1))
{ 
  print(i)
  z <- data.frame(CTRPv2auc1[,i],Yard[,"YardAUC"]) 
  z <- z[complete.cases(z),]
  z1 <- cor.test(z[,1],z[,2],method = "pearson")
  z1.cor[i] <- z1$estimate
  z1.pval[i] <- z1$p.value
}
names(z1.cor) <- colnames(CTRPv2auc1)
z1.cor1 <- data.frame(colnames(CTRPv2auc1),z1.cor)

final.CTRPv2auc <- z1.cor1[order(-z1.cor),]

pdf("RankingCTRPv2-Spearman.pdf")
plot(final.CTRPv2auc$z1.cor,col="red",xlab="Drug Index",ylab="Correlation",main="Correlation: Radiation and CTRPv2",cex.axis = 1.5,cex.lab = 2)
dev.off()

########################################################################

# PCL enrichment analysis using radiation AUCand SF2
# Step1: Correlation betweeen Radiation AUC and Chemo AUC across histologies
# Classes as gene sets and run enrichment analysis

source("plotRunningSum.R")

classes <- read.csv('CTRPv2Drugs_BasedOnGDSC1000.csv',stringsAsFactors=FALSE,sep="\t")
classes <- classes[,-c(5,6)]
classes$PossibleGDSC1000Class[grep(":",classes$PossibleGDSC1000Class)] <- "other"
classes$PossibleGDSC1000Class[which(classes$PossibleGDSC1000Class == "")] <- "other"
classes$PossibleGDSC1000Class[grep("chr",substring(classes$PossibleGDSC1000Class, 1, 3))] <- "chromatin"

genesets <- cbind(classes$Drugs,classes$PossibleGDSC1000Class)
colnames(genesets) <- c("Drugs","Class")

a <- loadGSC(genesets)

stats <- as.vector(final.CTRPv2auc$z1.cor)
names(stats) <- as.vector(final.CTRPv2auc$colnames.CTRPv2auc1.)
gsares <- runGSA(geneLevelStats=stats,gsc=a,nPerm=10000,geneSetStat="gsea",adjMethod="none")
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1AUC <- gseares1[,-(4:5)]
res <- data.frame(gseares1AUC$Name,gseares1AUC$`Stat (dist.dir)`,gseares1AUC$FDR)
res[,4] <- NA
res$V4[res$gseares1AUC.FDR<0.05] <- "red"
res$V4[res$gseares1AUC.FDR>0.05] <- "blue"

res1 <- res
colnames(res1) <- c("Name","ES","FDR","Color")
res1 <- res1[order(-res1$ES),]

res1[,5] <- NA
res1$V5[res1$FDR<0.05] <- "Significant(FDR<0.05)"
res1$V5[res1$FDR>0.05] <- "Non-significant"

library(ggplot2)
res1$Name <- factor(as.character(res1$Name), levels = as.character(res1$Name))
plt <- ggplot(res1, aes_string(x="Name", y= "ES", fill="V5") ) + geom_bar(stat = "identity")
  
plt <- plt + theme_bw() + theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))

##----remove x axis ------------------
plt <- plt +theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(), axis.line.x = element_blank())

##--------- add line at 0 ------------
plt <- plt + geom_hline(yintercept=0, size =0.25)
plt +  geom_text(aes(label=as.character(res1$Name)),size=4, angle = 90)
codY <- -1*res1$ES/abs(res1$ES)
hjustV <- ifelse((-1*codY)<0, 0,1)
codYV <- codY*0.05
plt <- plt + geom_text(y=codYV, label=as.character(res1$Name), hjust=hjustV, angle = 90)
plt

pdf("EnrichmentScore_PCL.pdf")
plt
dev.off()





