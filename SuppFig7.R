# Supplementary Figure 7

setwd("HypoxiaDrugs/")

#a <- PharmacoGx::availablePSets()
#CTRPv2 <- downloadPSet("CTRPv2")
#CTRPv2auc <- PharmacoGx::summarizeSensitivityProfiles(CTRPv2, sensitivity.measure='auc_published',verbose = TRUE)

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
  z1 <- cor.test(z[,1],z[,2],method = "spearman")
  z1.cor[i] <- z1$estimate
  z1.pval[i] <- z1$p.value
}
names(z1.cor) <- colnames(CTRPv2auc1)
z1.cor1 <- data.frame(colnames(CTRPv2auc1),z1.cor)

final.CTRPv2auc <- z1.cor1[order(-z1.cor),]

pdf("RankingCTRPv2-Spearman.pdf")
plot(final.CTRPv2auc$z1.cor,col="red",xlab="Drug Index",ylab="Correlation",main="Correlation: Radiation and CTRPv2",
     cex.axis = 1.5,cex.lab = 2)
dev.off()
