############################################################
#
# Ths code is used to fit the dose response data using the
# linear quadratic model implemented in the RadioGx package 
# Obtain best curves base don Rsquare
# Figure 2A-2C
#
############################################################

setwd("~/RadioGxAnalysis/")

require(gdata)
require(caTools)
require(RadioGx)
require(readxl)

# Yard et al. study
load('DoseResponse.RData')
D <- as.numeric(c("0","1","2","3","4","5","6","8","10"))
AUC <- as.numeric()
res <- list()

for(i in 1:nrow(dose.response))
{ 
  print(i)
  z <- unlist(dose.response[i,1:8])
  z <- c(1,z)
  res[[i]] <- linearQuadraticModel(D,z,family = "normal")
  AUC[i] <- computeAUC(D,z,area.type = "Fitted")
}
names(res) <- rownames(dose.response)

res1 <- matrix(NA,nrow(dose.response),5)
rownames(res1) <- rownames(dose.response)
colnames(res1) <- c("Yard-AUC","LQ-Alpha","LQ-Beta","LQ-AUC","LQ-Rsquare")
rownames(res1) <- rownames(dose.response)
res1[,"Yard-AUC"] <- dose.response$AUC
res1[,"LQ-AUC"] <- AUC
for(j in 1:length(res))
{
  res1[j,"LQ-Alpha"] <- unlist(res[[j]][1])
  res1[j,"LQ-Beta"] <- unlist(res[[j]][2])
  res1[j,"LQ-Rsquare"] <- attributes(res[[j]])[[which(names(attributes(res[[j]])) == "Rsquare")]]
}
res1[] <- sapply(res1, function(x) as.numeric(as.character(x))) # convert to numeric
save(res1,file="ClevelandLQFitVal.RData")

# Restrict cell lines that has > 0.6 RSquared value
res2 <- res1[which(res1[,5] >= 0.6),]
save(res2,file="ClevelandLQFitVal-GoodCurves.RData")

######################################################
# Figure 2A: Sample LQ model fitting using RadioGx 

sensitive <- "SK-ES-1" #Ewings_sarcoma-peripheral_primitive_neuroectodermal_tumour
resistant <- "SNU-245" # BILIARY_TRACT resistant cel lline

load('DoseResponse.RData')
load('ClevelandLQFitVal-GoodCurves.RData')
# plot resistant cell line
res1 <- res2
D <- as.numeric(c("0","1","2","3","4","5","6","8","10"))

LQ_resistant <- c(1,as.numeric(dose.response[which(resistant == rownames(dose.response)),1:8])) # plot resistant cell line
LQ_sensitive <- c(1,as.numeric(dose.response[which(sensitive == rownames(dose.response)),1:8])) # plot sensitive cell line

pdf("Example_LQFit1.pdf")
RadioGx:::doseResponseCurve(Ds=list("SNU-245" = D,"SK-ES-1" = D),
                            SFs=list("Resistant Cell Line" = LQ_resistant,"Sensitive Cell Line" = LQ_sensitive), 
                            plot.type="Both",legends.label = "rsquared",cex = 1.55,cex.main = 1.75,lwd = 2)
dev.off()

######################################################
# Figure 2B: Histogram of fitted AUC values

pdf("Hist-AUC.pdf")
hist(res2[,4],breaks=50,col = "red",xlab="AUC",cex.lab = 2,cex.axis=2,main="LQ Model AUC values")
dev.off()

################################################################
# Figure 2C: Compare clonogenic and proliferative assays from the same lab

Yard.Prol.Clon.auc <- read_excel('yard 2016 clonogenic proliferative.xlsx',1)
Yard.Prol.Clon.auc <- data.frame(Yard.Prol.Clon.auc)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"
a2 <- gsub(badchars, "", Yard.Prol.Clon.auc$Cell.Line)
Yard.Prol.Clon.auc$Cell.line <- a2

Yard.Prol.Clon.sf2 <- read_excel('yard 2016 clonogenic proliferative.xlsx',2)
Yard.Prol.Clon.sf2 <- data.frame(Yard.Prol.Clon.sf2)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf2$Cell.Line)
Yard.Prol.Clon.sf2$Cell.line <- a2

Yard.Prol.Clon.sf4 <- read_excel('yard 2016 clonogenic proliferative.xlsx',3)
Yard.Prol.Clon.sf4 <- data.frame(Yard.Prol.Clon.sf4)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf4$Cell.Line)
Yard.Prol.Clon.sf4$Cell.line <- a2

Yard.Prol.Clon.sf6 <- read_excel('yard 2016 clonogenic proliferative.xlsx',4)
Yard.Prol.Clon.sf6 <- data.frame(Yard.Prol.Clon.sf6)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf6$Cell.Line)
Yard.Prol.Clon.sf6$Cell.line <- a2

Yard.Prol.Clon.sf8 <- read_excel('yard 2016 clonogenic proliferative.xlsx',5)
Yard.Prol.Clon.sf8 <- data.frame(Yard.Prol.Clon.sf8)
a2 <- gsub(badchars, "", Yard.Prol.Clon.sf8$Cell.Line)
Yard.Prol.Clon.sf8$Cell.line <- a2

# Correlation
sf2 <- cor.test(Yard.Prol.Clon.sf2$Clonogenic,Yard.Prol.Clon.sf2$Proliferative)
sf4 <- cor.test(Yard.Prol.Clon.sf4$Clonogenic,Yard.Prol.Clon.sf4$Proliferative)
sf6 <- cor.test(Yard.Prol.Clon.sf6$Clonogenic,Yard.Prol.Clon.sf6$Proliferative)
sf8 <- cor.test(Yard.Prol.Clon.sf8$Clonogenic,Yard.Prol.Clon.sf8$Proliferative)
auc <- cor.test(Yard.Prol.Clon.auc$Clonogenic,Yard.Prol.Clon.auc$Proliferative)

t4 = as.numeric(c(sf2$estimate,sf4$estimate,sf6$estimate,sf8$estimate,auc$estimate))

# With error bars to the bar plot
F <- as.numeric(c(sf2$estimate,sf4$estimate,sf6$estimate,sf8$estimate,auc$estimate))
L <- as.numeric(c(min(sf2$conf.int),min(sf4$conf.int),min(sf6$conf.int),min(sf8$conf.int),min(auc$conf.int)))
U <- as.numeric(c(max(sf2$conf.int),max(sf4$conf.int),max(sf6$conf.int),max(sf8$conf.int),max(auc$conf.int)))

pdf("Corplot-YardProlClon.pdf",width=8,height=8)
par(oma=c(2,2,2,2)) # bottom,left,top,right
b <- barplot(F, main=" Yard Proliferative vs. Yard Clonogenic",beside=TRUE,
             names.arg=c("2 Gy","4 Gy","6 Gy","8 Gy","AUC"),las=0.2,ylab="Correlation",
             col=c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854"),ylim=c(0,1),cex.axis = 1.75,cex.names = 1.5,cex.lab = 1.75)
arrows(x0=b, x1=b, y0=L, y1=U, code=3, angle=90) 
dev.off()


# scatter plots
pdf("ScatterPlot-SF2.pdf")
plot(Yard.Prol.Clon.sf2$Clonogenic,Yard.Prol.Clon.sf2$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF2")
dev.off()

pdf("ScatterPlot-SF4.pdf")
plot(Yard.Prol.Clon.sf4$Clonogenic,Yard.Prol.Clon.sf4$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF4")
dev.off()

pdf("ScatterPlot-SF6.pdf")
plot(Yard.Prol.Clon.sf6$Clonogenic,Yard.Prol.Clon.sf6$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF6")
dev.off()

pdf("ScatterPlot-SF8.pdf")
plot(Yard.Prol.Clon.sf8$Clonogenic,Yard.Prol.Clon.sf8$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="SF8")
dev.off()

pdf("ScatterPlot-AUC.pdf")
plot(Yard.Prol.Clon.auc$Clonogenic,Yard.Prol.Clon.auc$Proliferative,col = "red",xlab="Clonogenic Assay",ylab="Proliferative Assay",main="AUC")
dev.off()











