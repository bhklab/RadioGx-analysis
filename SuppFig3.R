# Supplementary Figure 3

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

pdf("SF2-AUC-LTertiles.pdf")
plot(c1$SF2,c1$AUC,col="red",cex.lab = 2,cex.axis=2,xlab="SF2",ylab="AUC",main="Sensitive")
legend("topleft", legend=sprintf("Correlation= %.2f (p<e-16)",e1$estimate),bty="n",cex=1.5,horiz = TRUE)
dev.off()

pdf("SF2-AUC-ITertiles.pdf")
plot(c2$SF2,c2$AUC,col="red",cex.lab = 2,cex.axis=2,xlab="SF2",ylab="AUC",main="Intermediate")
legend("topleft", legend=sprintf("Correlation= %.2f (p<e-15)",e2$estimate),bty="n",cex=1.5,horiz = TRUE)
dev.off()

pdf("SF2-AUC-UTertiles.pdf")
plot(c3$SF2,c3$AUC,col="red",cex.lab = 2,cex.axis=2,xlab="SF2",ylab="AUC",main="Resistant")
legend("topleft", legend=sprintf("Correlation= %.2f (p<e-16)",e3$estimate),bty="n",cex=1.5,horiz = TRUE)
dev.off()


