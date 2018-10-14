# Supplementary Figure 6

a <- read_excel("Tissue Spec with Categories and totals.xlsx",1)
a <- data.frame(a)

pdf("PosEnriched.pdf")
hist(a$number.positive,col="red",xlab="Positively enriched pathways",main="",cex.lab=1.5,cex.axis=1.5,xlim=c(0,7))
dev.off()

pdf("NegEnriched.pdf")
hist(a$Number.negative,col="red",xlab="Negatively enriched pathways",main="",cex.lab=1.5,cex.axis=1.5)
dev.off()
