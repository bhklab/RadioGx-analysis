# Supplementary Figures
                         
# Supplementary Figure 2A, 2B

load('../data/ClevelandLQFitVal.RData')
pdf("../results/SuppFig2A-Hist-RSquared.pdf")
hist(res1[,5],breaks=50,col = "red",xlab="AUC",cex.lab = 2,cex.axis=2,main="LQ Model AUC values")
abline(v=0.6, col="black")
dev.off()

load('../data/ClevelandLQFitVal.RData')
load('../data/YardHistology.RData')

tmp <- YardHistology$Primarysite[match(rownames(res1),rownames(YardHistology))]
final_res1 <- data.frame(as.numeric(res1[,'LQ-Rsquare']),tmp)
colnames(final_res1) <- c("LQ-Rsquare","histlo")
rownames(final_res1) <- rownames(res1)

tt1 <- table(final_res1$histlo)
final_res2 <- subset(final_res1, final_res1[,2] %in% names(tt1[tt1 >= 10])) # minimum 10 cell lines per tissue
final_res2 <- final_res2[-which(final_res2[,2] == ""),]

df4 <- final_res2
long <- melt(df4, id.vars = "histlo")
bymedian <- with(df4, reorder(df4[,2], as.numeric(as.character(df4[,1])), median)) # order by median of oxic
long$histlo <- factor(bymedian)
labels1 <- c(levels(long$histlo))
pdf("../results/Supp2B-RSquare_Histology.pdf",width=10,height=10)
par(oma=c(15,10,5,5)) # bottom,left,top,right
boxplot(value ~ variable + histlo, data = long,las=2,col=rep(c("#66C2A5"),length(unique(df4$histlo))),
        at = seq(1,16),xaxt='n',ylab="RSquare (Goodness of fit)",
        cex.names = 1.75,cex.axis = 1.75,cex.lab=1.75)
axis(1,at = seq(1,16),labels = labels1,las=2,cex.axis = 1.75)
dev.off()

# Supplementary Figure 3
Yard.Prol.Clon.auc <- read_excel('/data/yard2016clonogenicproliferative.xlsx',1)
Yard.Prol.Clon.auc <- data.frame(Yard.Prol.Clon.auc)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"
a2 <- gsub(badchars, "", Yard.Prol.Clon.auc$Cell.Line)
Yard.Prol.Clon.auc$Cell.line <- a2

Yard.Prol.Clon.sf2 <- read_excel('/data/yard2016clonogenicproliferative.xlsx',2)
Yard.Prol.Clon.sf2 <- data.frame(Yard.Prol.Clon.sf2)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf2$Cell.Line)
Yard.Prol.Clon.sf2$Cell.line <- a2

Yard.Prol.Clon.sf4 <- read_excel('/data/yard2016clonogenicproliferative.xlsx',3)
Yard.Prol.Clon.sf4 <- data.frame(Yard.Prol.Clon.sf4)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf4$Cell.Line)
Yard.Prol.Clon.sf4$Cell.line <- a2

Yard.Prol.Clon.sf6 <- read_excel('/data/yard2016clonogenicproliferative.xlsx',4)
Yard.Prol.Clon.sf6 <- data.frame(Yard.Prol.Clon.sf6)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf6$Cell.Line)
Yard.Prol.Clon.sf6$Cell.line <- a2

Yard.Prol.Clon.sf8 <- read_excel('/data/yard2016clonogenicproliferative.xlsx',5)
Yard.Prol.Clon.sf8 <- data.frame(Yard.Prol.Clon.sf8)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf8$Cell.Line)
Yard.Prol.Clon.sf8$Cell.line <- a2

# scatter plots
pdf("../results/SuppFig3A-ScatterPlot-SF2.pdf")
plot(Yard.Prol.Clon.sf2$Clonogenic,Yard.Prol.Clon.sf2$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF2")
dev.off()

pdf("../results/SuppFig3B-ScatterPlot-SF4.pdf")
plot(Yard.Prol.Clon.sf4$Clonogenic,Yard.Prol.Clon.sf4$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF4")
dev.off()

pdf("../results/SuppFig3C-ScatterPlot-SF6.pdf")
plot(Yard.Prol.Clon.sf6$Clonogenic,Yard.Prol.Clon.sf6$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF6")
dev.off()

pdf("../results/SuppFig3D-ScatterPlot-SF8.pdf")
plot(Yard.Prol.Clon.sf8$Clonogenic,Yard.Prol.Clon.sf8$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF8")
dev.off()

pdf("../results/SuppFig3E-ScatterPlot-AUC.pdf")
plot(Yard.Prol.Clon.auc$Clonogenic,Yard.Prol.Clon.auc$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="AUC")
dev.off()

# Supplementary Figure 4

load('../data/GSEA-AUC.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
res.auc <- subset(gseares1,gseares1[,"FDR"] < 0.05) 
auc.gsea <- gsares

load('../data/GSEA-SF2.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]
gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(10000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
res.sf2 <- subset(gseares1,gseares1[,"FDR"] < 0.05) 
sf2.gsea <- gsares

rm(gsares)

common <- intersect(res.auc$Name,res.sf2$Name)
auc.common <- res.auc[match(common,res.auc$Name),]
sf2.common <- res.sf2[match(common,res.sf2$Name),]

#write.csv(common,file="CommonPath-May19-2019.csv")

#### Pie chart of EXTRA pathways that are enriched using AUC and not with SF2
diff_pathways <- read.csv('../data/CommonPath-May19-2019.csv',1)
diff_pathways <- diff_pathways[,c(2,3)]
diff_pathways <- diff_pathways[order(diff_pathways$Pathway),]
a <- table(diff_pathways$Pathway)
a <- cbind(names(a),a)

lbls<-rownames(a)
pct<-round(as.numeric(a[,2])/sum(as.numeric(a[,2]))*100)
lbls<-paste(lbls, pct)
lbls<-paste(lbls, "%", sep = "")

pdf("../results/SuppFig4-LeftPanel-CommonPathways-AUCSF2.pdf",width=20,height=10)
pie(as.numeric(a[,2]), labels = lbls,cex=1.5,col= c("#8DD3C7","#FFFFB3" ,"#FFED6F","#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462","red"), 
    main="Common Pathways (35) enriched by SF2, AUC",cex.main=2)
dev.off()

#### Pie chart of EXTRA pathways that are enriched using AUC and not with SF2
z <- setdiff(res.auc$Name,res.sf2$Name)
#write.csv(z,file="PathwaySpecificAUC-May19-2019.csv")

common_pathways <- read.csv('../data/PathwaySpecificAUC-May19-2019.csv')
common_pathways <- common_pathways[,c(2,3)]
common_pathways <- common_pathways[order(common_pathways$Pathway),]
a1 <- table(common_pathways$Pathway)
a1 <- cbind(names(a1),a1)

lbls<-rownames(a1)
pct<-round(as.numeric(a1[,2])/sum(as.numeric(a1[,2]))*100)
lbls<-paste(lbls, pct)
lbls<-paste(lbls, "%", sep = "")

pdf("../results/SuppFig4_RightPanel-PathwaysSpecific-AUC.pdf",width=20,height=10)
pie(as.numeric(a1[,2]), labels = lbls,cex=1.5,col=c("green","#8DD3C7","#FFFFB3","#FFED6F","#BEBADA","#80B1D3","red"), 
    main="Pathways (42) enriched specific to AUC",cex.main=2)
dev.off()

#### Supplementary Figure 6

load('../data/Oxic_Hypoxia5LQres.RData')
load('../data/CCLERNAseq.RData')

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
  tmp <- cor.test(ccle.rna1[i,],as.numeric(as.character(Yard1[,'Oxic'])),method="pearson")
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
  tmp <- cor.test(ccle.rna1[i,],as.numeric(as.character(Yard1[,'Hyp5'])),method="pearson")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(ccle.rna1)
metric1 <- data.frame(Gene,cor.val,p.val)
metric1$Gene <- as.character(metric1$Gene)
colnames(metric1) <- c("Gene","Correlation","Pvalue")
metric1[,4] <- p.adjust(metric1$Pvalue,method = "fdr",n = length(metric1$Pvalue))
metric1.hyp <- subset(metric1,metric1$V4<0.05)


res <- data.frame(metric$Gene,metric$Correlation,metric1$Correlation)
colnames(res) <- c("gene","Oxic","Hypoxic")
res[which(res[,2] < 0 & res[,3] < 0) ,4] <- "red"
res[which(res[,2] > 0 & res[,3] > 0) ,4] <- "red"
res[which(res[,2] < 0 & res[,3] > 0) ,4] <- "blue"
res[which(res[,2] > 0 & res[,3] < 0) ,4] <- "blue"

pdf("../results/SuppFig6-OxicHypoxic-FDR.pdf")
plot(res[,2],res[,3],col=res[,4],xlab="Oxic condition",ylab="Hypoxic condition",
     cex.lab=1.5,cex.axis=1.5)
legend("topleft", c("Genes that swapped directionality","Otherwise"), 
       xpd = TRUE, bty = "n", pch = c(0,0), col=c("blue","red"), cex = 1.2)

dev.off()



# Supplementary Figure 7 
a <- read_excel("../data/TissueSpecwithCategoriesandtotals.xlsx",1)
a <- data.frame(a)

pdf("../results/SuppFig7A-PosEnriched.pdf")
hist(a$number.positive,col="red",xlab="Positively enriched pathways",main="",cex.lab=1.5,cex.axis=1.5,xlim=c(0,7))
dev.off()

pdf("../results/SuppFig7B-NegEnriched.pdf")
hist(a$Number.negative,col="red",xlab="Negatively enriched pathways",main="",cex.lab=1.5,cex.axis=1.5)
dev.off()

# Supplementary Figure 8
load('/data/CTRPv2AUCPublishedVals.RData')
load('/data/DoseResponse.RData')
load('/data/ClevelandLQFitVal-GoodCurves.RData')

dose.response <- dose.response[rownames(res2),]

common <- intersect(rownames(dose.response),colnames(CTRPv2auc))
Yard <- dose.response[common,]
CTRPv2auc1 <- t(CTRPv2auc[,common])

z1.cor <- NULL
z1.pval <- NULL

for(i in 1:ncol(CTRPv2auc1))
{ 
  #print(i)
  z <- data.frame(CTRPv2auc1[,i],Yard[,"YardAUC"]) 
  z <- z[complete.cases(z),]
  z1 <- cor.test(z[,1],z[,2],method = "spearman")
  z1.cor[i] <- z1$estimate
  z1.pval[i] <- z1$p.value
}
names(z1.cor) <- colnames(CTRPv2auc1)
z1.cor1 <- data.frame(colnames(CTRPv2auc1),z1.cor)

final.CTRPv2auc <- z1.cor1[order(z1.cor),]

pdf("../results/SuppFig8-RankingCTRPv2-Spearman.pdf")
plot(final.CTRPv2auc$z1.cor,col="red",xlab="Drug Index",ylab="Correlation (Spearman)",main="Correlation: Radiation and CTRPv2",cex.axis = 1.5,cex.lab = 2)
arrows(545, 0.573, x1 = 520, y1 = 0.573, length = 0.25, angle = 30)
text(480,0.57,"VX-680")
arrows(544, 0.420, x1 = 520, y1 = 0.420, length = 0.25, angle = 30)
text(480,0.42,"Tenopiside")
dev.off()

print("Supplementary Figures Completed Successfully !")


                 