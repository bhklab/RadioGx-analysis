setwd("~/Documents/R Work For Project/Mutation Data DNA Repair/Data")

library(RadioGx)
library(tidyr)
library(plyr)
library(metaviz)
library(effsize)
library(tidyverse)
library(ggplot2)
cols=c("#66C2A5","#FC8D62","#8DA0CB") 
####Download mutation data from RadioGx

###Insert Script here 

##Export CSV files to run though VEP online tool 

##re-import those files once VEP categories added (By hand!! Sorry!)
file_list <- list.files(path="~/Documents/R Work For Project/Mutation Data DNA Repair/Data")
##B/c set wd, can use this as a file path also

d<-lapply(file_list, read.csv)
lapply(d, names) ##Now have all gene names data sets in one list
VEP_HIGH<-colnames(d[[]][24])##Correct spacing in manually entered column name
#str(d[[1]]) 


################################################################################################
################################################################################################
################################################################################################
#List of Yard cell lines with senesitivity data:
#Yard<-load("~/Documents/R Work For Project/Cleveland2.RData")
load("~/Documents/R Work For Project/ClevelandLQFitVal-GoodCurves.RData")

#myraw <- Cleveland@sensitivity$raw
res2<-as.data.frame(res2)
cell.lines<-toupper(rownames(res2))
#radauc <- summarizeSensitivityProfiles(Cleveland, sensitivity.measure = "AUC_recomputed")
#cell.lines<-toupper(colnames(radauc))
cell.lines<-gsub("-", "", cell.lines)
cell.lines<-gsub("/", "", cell.lines)
cell.lines<-gsub(" ", "", cell.lines)
cell.lines<-gsub("\\.", "", cell.lines)
cell.lines<-gsub("\\,", "", cell.lines)
cell.lines<-gsub("\\(", "", cell.lines)
cell.lines<-gsub(")", "", cell.lines)

AUC<-as.vector(res2$`LQ-AUC`)

##in cell line names within the mutation data, get rid of the histology after_
for (i in 1:length(d)){
  d.n<-gsub("_.*", "", d[[i]][,3])  ##Have made the list of names we want to replace them with
  d[[i]][3]<-d.n
}
d.cl<-lapply(d, subset, Tumor.Sample.Barcode %in% cell.lines)

d.mut.high<-lapply(d.cl, subset, VEP_HIGH==TRUE)
d.mut.notHIGH<-lapply(d.cl, subset, VEP_HIGH==FALSE)

##Make sure none of the cell lines are both high and low, also eliminate duplicates
for(i in 1:length(d.mut.notHIGH)){
  n<-which(d.mut.notHIGH[[i]]$Tumor.Sample.Barcode %in% d.mut.high[[i]]$Tumor.Sample.Barcode)
  #print(n) ##uncomment if you want to see how many changes are made
  if(length(n)>0){
  d.mut.notHIGH[[i]]<-d.mut.notHIGH[[i]][-n,]
  }
  d.mut.notHIGH[[i]]<-distinct(d.mut.notHIGH[[i]], Tumor.Sample.Barcode, .keep_all = TRUE)
  #print(nrow(d.mut.notHIGH[[i]]))
}

##remove duplicates from High category too
for(i in 1:length(d.mut.high)){
  d.mut.high[[i]]<-distinct(d.mut.high[[i]], Tumor.Sample.Barcode, .keep_all = TRUE)
  print(nrow(d.mut.high[[i]]))
}


#Get the names of the genes
gene.names<-(gsub("_.*", "", file_list))

##now set up a data frame that incorporates all of the information 
##genes are columns and cell lines are rows 

##Make the matrix we will make boxplots from
Final<-matrix(data = 0, ncol=14, nrow=length(cell.lines), dimnames = list(cell.lines, gene.names))
###Make all the TRUE values in mut_notHIGH as "Low" in final matrix
L<-c()
for(i in 1:length(d.mut.notHIGH)){
  
  for(j in 1:nrow(d.mut.notHIGH[[i]])){
    
    k<-which(rownames(Final)== d.mut.notHIGH[[i]][j,3])
    L<-c(L,k)
    Final[k,i]<-"Low"
  }
  M<-length(L)
  print(paste(i, M, sep = " "))
  L<-c()
}

###Make all the TRUE values in mut_HIGH as 2 in final matrix 
L<-c()
for(i in 1:length(d.mut.high)){
  
  for(j in 1:nrow(d.mut.high[[i]])){
    
    k<-which(rownames(Final)== d.mut.high[[i]][j,3])
    L<-c(L,k)
    Final[k,i]<-"High"
  }
  M<-length(L)
  print(paste(i, M, sep = " "))
  L<-c()
  
}
Final<-replace(Final, Final==0, "WT")
Final<-as.data.frame(Final)
Final$Sensitivity<-AUC

##############################################################
##############################################################
##Repeat what we just did, but including information that incorporates 
#how many samples there are in each category 
Final.Labels<-matrix(data = 0, ncol=14, nrow=length(cell.lines), dimnames = list(cell.lines, gene.names))

L<-c()
B<-list() #This will be the list of which rows in the matrix are
##in the selected category
for(i in 1:length(d.mut.notHIGH)){
  
  for(j in 1:nrow(d.mut.notHIGH[[i]])){
    
    k<-which(rownames(Final.Labels)== d.mut.notHIGH[[i]][j,3])
    L<-c(L,k)
  }
  Final.Labels[L,i]<-paste("VEP Not High, n=",length(L), sep = "" )
  print(paste(i, length(L), sep = " ")) ##Show how many mutations/gene... just for interest
  L<-list(L)
  B[i]<-L
  L<-c()
}

###Make all the TRUE values in mut_HIGH as "HIGH" in final matrix 
L<-c()
P<-list()#This will be the list of which rows in the matrix are
##in the selected category
for(i in 1:length(d.mut.high)){
  for(j in 1:nrow(d.mut.high[[i]])){
    k<-which(rownames(Final.Labels)== d.mut.high[[i]][j,3])
    L<-c(L,k)
  }
  Final.Labels[L,i]<-paste("VEP High, n=",length(L), sep = "" )
  print(paste(i, length(L), sep = " "))
  L<-list(L)
  P[i]<-L
  L<-c()
}
#Get the same list/change names for WT 
L<-c()
K<-list()
for(i in 1:ncol(Final.Labels)) {
  k<-as.numeric(which(Final.Labels[,i]== 0))
  L<-c(L,k)
  Final.Labels[L,i]<-paste("WT, n=",length(L), sep = "" )
  print(paste(i, length(L), sep = " "))
  L<-list(L)
  K[i]<-L
  L<-c()
}
Final.Labels<-as.data.frame(Final.Labels)
Final.Labels$Sensitivity<-AUC
Final.Labels2<- gather(Final.Labels, Gene, Status, #make this into a long format that ggplot2 understands
                       ATM:XRCC5, factor_key=TRUE) 
x.labels<-Final.Labels2$Status
##############################################################################################

#calculate cohen d values of effect size for all of the mutations,
##And also calculate wilcoxon p value (FDR adjust) for each gene
cohen.values<-c()
cohen.se<-c()
cohen.rowname<-c()
name.p.values<-c()
for(i in 1:14){
  
  T<-Final.Labels[P[[i]],c(i,15)] 
  F<-Final.Labels[K[[i]],c(i,15)]
  O<-Final.Labels[B[[i]], c(i,15)]
  m<- T$Sensitivity      ###For cohen values - high mutated
  wt<- F$Sensitivity       ####For cohen values - wildtype
  
  d<-c(m,wt)
  n1<-length(m)
  n2<-length(wt)
  df<-n1+n2-1
  f<-rep(c("mutated", "wildtype"), c(n1, n2))
  cohen.temp<-cohen.d(d~f)
  test.sd<-sd(T$Sensitivity)
  test.se<-sqrt(((n1+n2)/(n1*n2) + cohen.temp$estimate/(2*df))*((n1+n2)/df)) #Formula to get variance of cohen values
  #temp.se<-qnorm(0.975)*test.sd/sqrt(length(T$Sensitivity))
  
  cohen.rowname<-c(cohen.rowname, paste(colnames(Final.Labels[i]), ", n=", length(P[[i]]), sep =""))
  cohen.values<-c(cohen.values, cohen.temp$estimate) ##testing - switched out from temp.se
  cohen.se<-c(cohen.se, test.se)
  w<-wilcox.test(T$Sensitivity,F$Sensitivity)
  p.value<-round(p.adjust(w$p.value, method = "bonferroni", n=14), digits = 6)
  name.p.values<-c(name.p.values, paste(colnames(Final.Labels[i]), ", p=", p.value, sep="")) #Adds to list of p values

}
rm(P)
rm(K)
rm(B)

##############################################################################################
##############################################################################################

colnames(Final)<-c(name.p.values, "Sensitivity")

Final2<-gather(Final, Gene, Status, "ATM, p=0.114978":"XRCC5, p=1", factor_key=TRUE) #I manually made these names fit
Final2$x.labels<-x.labels

Individual.Genes<- ggplot(data = Final2, aes(x=x.labels, y=Sensitivity )) +
  geom_boxplot(aes(fill=factor(Status)), outlier.shape = NA)+facet_wrap( ~ Gene, scales="free")+
  scale_fill_manual(values = cols[1:3])+
  scale_y_continuous(name = "" )+ #removed AUC for omnigraffle image
  labs(x="")+
  theme_linedraw(base_size=26)+ ##Massive because there are so many boxes it ends up being tiny
  theme(legend.position = "none", axis.text.x = element_text(size= 12),
        axis.text.y = element_text(size= 18))+
  geom_jitter(data=Final2, size=1, shape = 1, position=position_jitter(width=.2)) 

Individual.Genes

##############################################################################################
##############################################################################################

##Making the total mutations plot 

Final$TotalHigh <- rowSums(Final == "High") ##This column refers to high mutations
Final$TotalLow<-rowSums(Final =="Low") ##Recall low just means not VEP High 

N<-which(Final$TotalHigh == 0) ##N = Not High 
M<-which(Final$TotalLow >0 & Final$TotalHigh ==0) ##Which only low 
Either<-which(Final$TotalHigh>0 | Final$TotalLow>0) ##Which has one or the other 

Final$AnyMutationStatus<-NA #Make a new column

##15 is the AnyMutationStatus column number
med.mut<-Final[M, 15] ##Medium mutations! (VEP not high)
yes.mut<-Final[-N, 15] ## recall N = not high, so -N is High 
no.mut<-Final[-Either, 15]

Final$AnyMutationStatus[-N]<-paste("VEP High, n=",length(yes.mut), sep = "" )
Final$AnyMutationStatus[M]<-paste("VEP Not High, n=", length(med.mut), sep = "")
Final$AnyMutationStatus[-Either]<-paste("WT, n=", length(no.mut), sep = "" )
w<-wilcox.test(yes.mut, no.mut)
p.value<-round(w$p.value,digits=4)

total.plot<-ggplot(data = Final, aes(x = AnyMutationStatus, y = Sensitivity))+
  geom_boxplot(aes(fill=factor(AnyMutationStatus)), outlier.shape = NA)+ 
  labs(x="")+
  scale_y_continuous(name = "")+ ##removed - will add label in omnigraffle
  scale_fill_manual(values = cols[1:3], guide = FALSE)+
  theme_linedraw(base_size=45)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(paste("All Genes, p=", p.value, sep = ""))+
geom_jitter(data=Final, size=1, shape = 1, position=position_jitter(width=.2))
  

total.plot
##############################################################################################
##############################################################################################

##Forest Plot
##All these components were made earlier when making supp figure 


matrix.col<-c("Cohen", "SE")
Forest.matrix<-matrix(c(cohen.values, cohen.se),nrow = 14, ncol=2)
Forest.matrix<-as.data.frame(Forest.matrix)
rownames(Forest.matrix)<-cohen.rowname
colnames(Forest.matrix)<-matrix.col
as.matrix<-Forest.matrix
Forest.matrix$Gene<-rownames(Forest.matrix)

Forest.matrix<-Forest.matrix %>% arrange(Cohen)
rownames(Forest.matrix)<-Forest.matrix$Gene
Forest.matrix<-Forest.matrix[,-3]

viz_thickforest(Forest.matrix, study_labels = rownames(Forest.matrix), text_size= 5,
               col=cols[1],tick_col=cols[2], summary_col = cols[3],
                  xlab = "Cohen d", annotate_CI = FALSE) 

