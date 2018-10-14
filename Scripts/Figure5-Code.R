##################################################
#
# Normoxic = 21% of 760 Atmospheric Pressure
# 21% is Normoxia in in-vitro,0.65%
#
##################################################

# Pathway analysis carried out using YardAUC and modelled hypoxia LQAUC values

library(RadioGx)
library(xlsx)
library(gdata)
library(piano)

setwd("~/RadioGxAnalysis/")

##### Mathematical modelling of hypoxia
# Fig5A

OER_m = 3
K_m = 3
po2 = 160
a = ((OER_m*po2)+K_m)/(po2+K_m)
OMF = (1/OER_m)*a

D <- as.numeric(c("0","1","2","3","4","5","6","8","10"))
SF1 = exp(-0.3*D*OMF-(0.03*D*D*OMF))
RadioGx::computeAUC(D,SF1)

po2 = 10
a = ((OER_m*po2)+K_m)/(po2+K_m)
OMF = (1/OER_m)*a
SF2 = exp(-0.3*D*OMF-(0.03*D*D*OMF))
RadioGx::computeAUC(D,SF2)

po2 = 5
a = ((OER_m*po2)+K_m)/(po2+K_m)
OMF = (1/OER_m)*a
SF3 = exp(-0.3*D*OMF-(0.03*D*D*OMF))
RadioGx::computeAUC(D,SF3)

pdf("HyxpoxiaPlot.pdf")
RadioGx:::doseResponseCurve(Ds=list("Normoxia (21%=160mmHG)" = D,"Hypoxia (1.3%=10mmHG)" = D,"Hypoxia (0.65%=5mmHG)" = D),
                            SFs=list("Hypoxia" = SF1,"Hypoxia" = SF2,"Hypoxia" = SF3), plot.type="Actual",legends.label = NULL,title = "Effect of Hypoxia",
                            cex = 1.55,cex.main = 1.75,lwd = 2)
dev.off()

# Fig 5B
# load the LQ model values for good curves
load('DoseResponse.RData')
load('ClevelandLQFitVal-GoodCurves.RData')
dose.response <- res2

D <- as.numeric(c("0","1","2","3","4","5","6","8","10"))

# {Oxic,Hypoxic} = {60,5} and Fractionation = 1
hypoxia <- as.numeric(c("5"))
hypoxia.res <- list()
for(j in 1:length(hypoxia))
{ 
  print(j)
  KR = 1
  OER_m = 3
  K_m = 3
  po2 = hypoxia[j]
  a = ((OER_m*po2)+K_m)/(po2+K_m)
  OMF = (1/OER_m)*a
  SF1AUC <- as.numeric()
  for(i in 1:nrow(dose.response))
  { 
    print(i)
    SF1 <- exp(-dose.response[i,'LQ-Alpha']*D*OMF-((dose.response[i,'LQ-Beta']*D*D*OMF*OMF)/KR))
    SF1AUC[i] <- computeAUC(D,SF1)
  }
  hypoxia.res[[j]] <- SF1AUC
}
hypoxia.LQ <- do.call(cbind, hypoxia.res)
rownames(hypoxia.LQ) <- rownames(dose.response)
colnames(hypoxia.LQ) <- c("Hyp5")
test <- cbind(res2[,'Yard-AUC'],hypoxia.LQ) 
colnames(test)[1] <- c("Oxic")
hypoxia.LQ <- test
save(hypoxia.LQ,file="Oxic_Hypoxia5LQres.RData")

######################################################

# Pathway Enrichment Analysis: O2 = {YardAUC,5mmHg} = {21% is Normoxia in in-vitro,0.65%}
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
res.gsea <- subset(gseares1,gseares1[,"FDR"] < 0.1) 
#save(gsares,file="GSEA-Oxic.RData")
save(gsares,file="GSEA-Hyp5.RData")

##########################################################################################
# Oxic and Hypoxic box plots side by side

load('Oxic_Hypoxia5LQres.RData')
load('YardHistology.RData')

histlo <- YardHistology$Primarysite[match(rownames(hypoxia.LQ),rownames(YardHistology))]
test <- data.frame(hypoxia.LQ,histlo)
test$Oxic <- as.numeric(as.character(test$Oxic))
test$Hyp5 <- as.numeric(as.character(test$Hyp5))
test$histlo <- (as.character(test$histlo))

df <- data.frame(test)
df$histlo <- as.character(df$histlo)
df$histlo[which(df$histlo == "")] <- "Others"

tt <- table(df$histlo)
df2 <- subset(df, df$histlo %in% names(tt[tt > 15]))

df3 <- df2[-which(df2$histlo == "Others"),]

long <- melt(df3, id.vars = "histlo")
bymedian <- with(df3, reorder(df3$histlo, df3$Oxic, median)) # order by median of oxic
long$histlo <- factor(bymedian)

labels1 <- c("Endometrium","Stomach","Large Intestine","Ovary","Upper AeroTract","Lung","Oesophagus",
             "Pancreas","Urinary Tract","Liver","CNS","Skin","Breast")

pdf("CompareAUC-OxicHyp.pdf",width=10,height=10)
par(oma=c(15,10,5,5)) # bottom,left,top,right
boxplot(value ~ variable + histlo, data = long,las=2,col=rep(c("#66C2A5","#FC8D62"),length(unique(df3$histlo))),
        at = c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38),xaxt='n',ylab="AUC Value",
        cex.names = 1.75,cex.axis = 1.75)
axis(1,at = seq(1,39,by=3)+0.5,labels = labels1,las=2,cex.axis = 1.75)
legend("bottom", c("Oxic","Hypoxic"), xpd = TRUE, horiz = TRUE, inset = c(1.0,1), bty = "n", pch = c(0,0), col=c("#66C2A5","#FC8D62"), cex = 1.5)
dev.off()


##################################################################################

# Figure 5C

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

metric2 <- data.frame(metric$Gene,metric$Correlation-metric1$Correlation)

oxic <- metric[order(-metric$Correlation),]
oxic[,4] <- rep(1:nrow(oxic))
colnames(oxic)[4] <- "Rank"

hypoxic <- metric1[order(-metric1$Correlation),]
hypoxic[,4] <- rep(1:nrow(hypoxic))
colnames(hypoxic)[4] <- "Rank"

hypoxic1 <- hypoxic[match(oxic$Gene,hypoxic$Gene),]

final.AUC <- oxic
final.SF2 <- hypoxic
# change in rankings
rankAUC <- final.AUC 
rankAUC[,"Index"] <- rep(1:nrow(rankAUC))
rankSF2 <- final.SF2
rankSF2[,"Index"] <- rep(1:nrow(rankSF2))
deltarank <- matrix(NA,nrow=nrow(rankAUC),2)
deltarank[,1] <- as.character(rankAUC$Gene)
deltarank[,2] <- abs(rankAUC$Index - match(rankAUC$Gene,rankSF2$Gene))
colnames(deltarank) <- c("Name","RankChange")
rownames(deltarank) <- deltarank[,'Name']
deltarank <- data.frame(deltarank)
deltarank$RankChange <- as.numeric(as.character(deltarank$RankChange))

# plot delta rank change by index
test <- deltarank[order(deltarank$RankChange),]
pdf("DeltaChange-OxicHypoxic.pdf")
plot(test$RankChange,xlab="Index of Gene",ylab="Delta Rank",col="red")
dev.off()




