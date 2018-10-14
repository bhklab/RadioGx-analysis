source("plotRunningSum.R")

load('GSEA-AUC.RData')
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

load('GSEA-SF2.RData')
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

# extract leading edge genes
leadingedgegenes.auc <- list()
GeneSets.auc <- list()
for(k in 1:nrow(auc.common)){
  t <- as.character(auc.common$Name[k])
  par(mar = rep(2, 4))
  aa <- plotRunningSum(gsaRes=auc.gsea,geneSet=t,1) # Plot the enrichment score
  leadingedgegenes.auc[[k]] <- aa$leadingEdge
  GeneSets.auc[[k]] <- names(leadingedgegenes.auc[[k]])
}
names(GeneSets.auc) <- auc.common$Name

leadingedgegenes.sf2 <- list()
GeneSets.sf2 <- list()
for(k in 1:nrow(sf2.common)){
  t <- as.character(sf2.common$Name[k])
  par(mar = rep(2, 4))
  aa <- plotRunningSum(gsaRes=sf2.gsea,geneSet=t,1) # Plot the enrichment score
  leadingedgegenes.sf2[[k]] <- aa$leadingEdge
  GeneSets.sf2[[k]] <- names(leadingedgegenes.sf2[[k]])
}
names(GeneSets.sf2) <- sf2.common$Name

# intersection of genes
common.leg <- as.numeric()
for(i in 1:length(GeneSets.auc)){
  common.leg[i] <- length(intersect(GeneSets.sf2[[i]],GeneSets.auc[[i]]))
}

leadingedgegenes <- data.frame(sapply(GeneSets.auc, length),sapply(GeneSets.sf2, length),common.leg)
colnames(leadingedgegenes) <- c("AUC","SF2","Common")
rownames(leadingedgegenes) <- sf2.common$Name
d <- do.call(rbind, leadingedgegenes)

pdf("LeadingEdgeGenes.pdf",width=20,height=10)
par(mar = c(32, 4, 4, 2) + 0.1)
barplot(d, beside = TRUE,args.legend = list(x = "topleft", bty="n"),col=c("gold3","indianred1","grey"),ylab="Number of Leading Edge Genes",las=2,
        names.arg= auc.common$Name,cex.axis=1, cex.names=1)
legend("topleft", c("AUC-LeadingEdgeGenes","SF2-LeadingEdgeGenes","Common"), pch=15, col=c("gold3","indianred1","grey"), bty="n")
dev.off()


#### Pie chart of EXTRA pathways that are enriched using AUC and not with SF2
z <- setdiff(res.auc$Name,res.sf2$Name)
diff_pathways <- read.csv('Pathwaycategories-Diff_AUCSF2 Update Jan14.csv',1)
diff_pathways <- diff_pathways[,c(2,3)]
diff_pathways <- diff_pathways[order(diff_pathways$Pathway),]
a <- table(diff_pathways$Pathway)
a <- cbind(names(a),a)

lbls<-rownames(a)
pct<-round(as.numeric(a[,2])/sum(as.numeric(a[,2]))*100)
lbls<-paste(lbls, pct)
lbls<-paste(lbls, "%", sep = "")

pdf("DiffPathways-AUCSF2.pdf",width=20,height=10)
pie(as.numeric(a[,2]), labels = lbls,cex=1.5,col= c("#8DD3C7","#FFFFB3" ,"#FFED6F","#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462"), 
    main="Pathways (44) enriched specific to AUC",cex.main=2)
dev.off()

#### Pie chart of EXTRA pathways that are enriched using AUC and not with SF2
#z <- intersect(res.auc$Name,res.sf2$Name)
common_pathways <- read.csv('Pathways Categories AUCSF2 Update Jan14.csv')
common_pathways <- common_pathways[,c(2,3)]
common_pathways <- common_pathways[order(common_pathways$Pathway),]
a1 <- table(common_pathways$Pathway)
a1 <- cbind(names(a1),a1)

lbls<-rownames(a1)
pct<-round(as.numeric(a1[,2])/sum(as.numeric(a1[,2]))*100)
lbls<-paste(lbls, pct)
lbls<-paste(lbls, "%", sep = "")

pdf("CommonPathways-AUCSF2.pdf",width=20,height=10)
pie(as.numeric(a[,2]), labels = lbls,cex=1.5,col=c("#8DD3C7","#FFFFB3" ,"#FFED6F","#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462" ,"#B3DE69"), main="Common Pathways (39) enriched using AUC, SF2",cex.main=2)
dev.off()
