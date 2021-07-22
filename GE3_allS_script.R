#1.1
#Creation of a matrix that includes the given data
gednorm<-as.matrix(read.delim("GeneExpressionDataset_normalized.tsv",header=T,row.names=1))

#Creation of a matrix where the logFC values and the p-values will be stored
genevalues<-matrix(NA,length(gednorm[,1]),10,dimnames=list(rownames(gednorm),c("logFC(TG)","p-value(TG)","logFC(TherA)","p-value(TherA)","logFC(TherB)","p-value(TherB)","logFC(TherC)","p-value(TherC)","logFC(TherD)","p-value(TherD)")))

#Performance of a for loop for every given gene that creates a dataframe in every run with the 
#proper settings on which the ANOVA analysis will take place. Then, the needed values are stored
#to the "genevalues" matrix for each gene
for(i in 1:length(gednorm[,1])){
  geneexp<-data.frame(matrix(NA,60,2))
  colnames(geneexp)<-c("Expression","Condition")
  geneexp[,2]<-c(rep("WT",10),rep("TG",10),rep("TherA",10),rep("TherB",10),rep("Therc",10),rep("TherD",10))
  geneexp[,1]<-unlist(unname(gednorm[i,]))
  fit<-aov(Expression~Condition,data=geneexp)
  results<-TukeyHSD(fit)
  genevalues[i,1]<- -(results[[1]][5])
  genevalues[i,2]<- (results[[1]][50])
  genevalues[i,3]<- -(results[[1]][9])
  genevalues[i,4]<- (results[[1]][54])
  genevalues[i,5]<- -(results[[1]][12])
  genevalues[i,6]<- (results[[1]][57])
  genevalues[i,7]<- -(results[[1]][14])
  genevalues[i,8]<- (results[[1]][59])
  genevalues[i,9]<- -(results[[1]][15])
  genevalues[i,10]<- (results[[1]][60])
}
head(genevalues)

#Calculation of the amount of differentially expressed genes and creation of a matrix with 
#the appopriate number of rows
x=0
y=1
for (i in 1:length(genevalues[,1])){
  if(((abs(genevalues[i,1])>=1)&(genevalues[i,2]<=0.05))|((abs(genevalues[i,3])>=1)&(genevalues[i,4]<=0.05))|((abs(genevalues[i,5])>=1)&(genevalues[i,6]<=0.05))|((abs(genevalues[i,7])>=1)&(genevalues[i,8]<=0.05))|((abs(genevalues[i,9])>=1)&(genevalues[i,10]<=0.05))){
    x<-x+1
  }}
DEG_list<-matrix(NA,x,1)
#Performance of a for and an if loop that check which genes are differentially expressed
#(according to the given criteria) and store their names to the "DEG_list" matrix
for (i in 1:length(genevalues[,1])){
  if(((abs(genevalues[i,1])>=1)&(genevalues[i,2]<=0.05))|((abs(genevalues[i,3])>=1)&(genevalues[i,4]<=0.05))|((abs(genevalues[i,5])>=1)&(genevalues[i,6]<=0.05))|((abs(genevalues[i,7])>=1)&(genevalues[i,8]<=0.05))|((abs(genevalues[i,9])>=1)&(genevalues[i,10]<=0.05))){
    DEG_list[y,1]<-rownames(genevalues)[i]
    y<-y+1
  }}
library("openxlsx")
write.xlsx(DEG_list,"1.DEG_list.xlsx",col.names=F)

#1.2
#Creation of a matrix that contains the positions that hold the "DEG_list" genes in the "genevalues" matrix
z<-matrix(NA,length(DEG_list),1)
for(i in 1:length(DEG_list)){
  z[i,1]<-which(rownames(genevalues)==DEG_list[i,1])}
head(z)

#Creation of a matrix that contains only the logFC values of the differentially expressed genes
heatdata<-matrix(NA,length(DEG_list),5)
rownames(heatdata)<-DEG_list
colnames(heatdata)<-c("TG","TherA","TherB","TherC","TherD")
for (i in 1:length(z)){
  heatdata[i,1]<-genevalues[z[i,1],1]
  heatdata[i,2]<-genevalues[z[i,1],3]
  heatdata[i,3]<-genevalues[z[i,1],5]
  heatdata[i,4]<-genevalues[z[i,1],7]
  heatdata[i,5]<-genevalues[z[i,1],9]
}
head(heatdata)

#Activation of the needed libraries to visualize the heatmap with a bit more interesting colours
library(gplots)
library("RColorBrewer")
cols<- colorRampPalette(brewer.pal(11, "Spectral"))(25)
heatmap.2(heatdata,col=cols,cexCol=1.2,main="Differentially Expressed Genes")

#1.3
#Use of the "NbClust" function that performs tests and ends up with the best number of k clusters
#for the given data
library(NbClust)
nb<-NbClust(heatdata,min.nc=5,max.nc=8,method="ward.D2")

#Performance of the "kmeans" function seperating the genes in the suggested by the "NbClust"
#amount of clusters
set.seed(93)
kclusters<-kmeans(heatdata,centers=5,iter.max=100)
cluster1<-which(kclusters$cluster==1)
cluster2<-which(kclusters$cluster==2)
cluster3<-which(kclusters$cluster==3)
cluster4<-which(kclusters$cluster==4)
cluster5<-which(kclusters$cluster==5)

#Creation of the character vectors that contain the genes of each clusters
cluster1.genenames<-rownames(as.matrix(cluster1))
cluster2.genenames<-rownames(as.matrix(cluster2))
cluster3.genenames<-rownames(as.matrix(cluster3))
cluster4.genenames<-rownames(as.matrix(cluster4))
cluster5.genenames<-rownames(as.matrix(cluster5))
write.xlsx(cluster1.genenames,"3.Genenames from cluster 1.xlsx",col.names=F)
write.xlsx(cluster2.genenames,"3.Genenames from cluster 2.xlsx",col.names=F)
write.xlsx(cluster3.genenames,"3.Genenames from cluster 3.xlsx",col.names=F)
write.xlsx(cluster4.genenames,"3.Genenames from cluster 4.xlsx",col.names=F)
write.xlsx(cluster5.genenames,"3.Genenames from cluster 5.xlsx",col.names=F)

# 1.4
#Installation and activation of the "gProfileR" package that downloads functional analysis data
#from given databases and attaches them to the appropriate genes
install.packages("gProfileR", repos = "http://cran.us.r-project.org")
library(gProfileR)
funcenr1<-gprofiler(query=as.character(cluster1.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr2<-gprofiler(query=as.character(cluster2.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr3<-gprofiler(query=as.character(cluster3.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr4<-gprofiler(query=as.character(cluster4.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr5<-gprofiler(query=as.character(cluster5.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))

#Creation of character vectors that contain the "term.name" and "p.value" values of the 
#most significant genes according to given criteria
cluster1.term.name<-funcenr1$term.name[which((funcenr1$term.size<=200)&(funcenr1$p.value<=0.01))]
cluster2.term.name<-funcenr2$term.name[which((funcenr2$term.size<=200)&(funcenr2$p.value<=0.01))]
cluster3.term.name<-funcenr3$term.name[which((funcenr3$term.size<=200)&(funcenr3$p.value<=0.01))]
cluster4.term.name<-funcenr4$term.name[which((funcenr4$term.size<=200)&(funcenr4$p.value<=0.01))]
cluster5.term.name<-funcenr5$term.name[which((funcenr5$term.size<=200)&(funcenr5$p.value<=0.01))]

cluster1.pvalue<-funcenr1$p.value[which((funcenr1$term.size<=200)&(funcenr1$p.value<=0.01))]
cluster2.pvalue<-funcenr2$p.value[which((funcenr2$term.size<=200)&(funcenr2$p.value<=0.01))]
cluster3.pvalue<-funcenr3$p.value[which((funcenr3$term.size<=200)&(funcenr3$p.value<=0.01))]
cluster4.pvalue<-funcenr4$p.value[which((funcenr4$term.size<=200)&(funcenr4$p.value<=0.01))]
cluster5.pvalue<-funcenr5$p.value[which((funcenr5$term.size<=200)&(funcenr5$p.value<=0.01))]

#Creation of the dataframes that contain the gene
ssf1<-data.frame(cluster1.term.name,cluster1.pvalue)
ssf2<-data.frame(cluster2.term.name,cluster2.pvalue)
ssf3<-data.frame(cluster3.term.name,cluster3.pvalue)
ssf4<-data.frame(cluster4.term.name,cluster4.pvalue)
ssf5<-data.frame(cluster5.term.name,cluster5.pvalue)

#Arranging the dataframe with increasing order with the "arrange" function of the "tidyverse" library
library(tidyverse)
ssf1<-ssf1%>%arrange(cluster1.pvalue)
ssf2<-ssf2%>%arrange(cluster2.pvalue)
ssf3<-ssf3%>%arrange(cluster3.pvalue)
ssf4<-ssf4%>%arrange(cluster4.pvalue)
ssf5<-ssf5%>%arrange(cluster5.pvalue)
write.xlsx(ssf1,"4.Cluster 1 Functions.xlsx",col.names=T)
write.xlsx(ssf2,"4.Cluster 2 Functions.xlsx",col.names=T)
write.xlsx(ssf3,"4.Cluster 3 Functions.xlsx",col.names=T)
write.xlsx(ssf4,"4.Cluster 4 Functions.xlsx",col.names=T)
write.xlsx(ssf5,"4.Cluster 5 Functions.xlsx",col.names=T)

#1.5
#Creation of a matrix and manipulation of the data so as to be inserted to the "randomForest" function
tRFdata<-matrix(NA,60,length(DEG_list[,1]))
colnames(tRFdata)<-DEG_list
State<-factor(c(rep("WT",10),rep("TG",10),rep("TherA",10),rep("TherB",10),rep("TherC",10),rep("TherD",10)),levels=c("WT","TG","TherA","TherB","TherC","TherD"))
for(i in 1:length(tRFdata[1,])){
  tRFdata[,i]<-gednorm[(which(rownames(gednorm)==DEG_list[i,1])),]}
tRFdata<-as.data.frame(cbind.data.frame(tRFdata,State))
View(tRFdata)

#Activation of the "randomForest" library
library(randomForest)
#Conversion of all the names to valid due to some errors that kept appearing due to 
#some symbols that they contained
names(tRFdata)<-make.names((names(tRFdata)))
set.seed(85)
RFmodel<-randomForest(State~.,data=tRFdata,ntree=1000,mtry=3)
ConfusionMatrix<-RFmodel$confusion
write.xlsx(ConfusionMatrix,"5.Confusion Matrix.xlsx",col.names=T,row.names=T)


#1.6
#Creation of a matrix where the data from the "varImp" function will be stored
library(caret)
varib<-as.matrix(varImp(RFmodel))
#Arranging of the "varib.ord" matrix to decreasing order so as to seperate the 50 most important genes
varib.ord<-as.matrix(varib[order(varib,decreasing=T),])
#Seperation of the 50 most important genes to the "imp.genes" matrix
imp.genes<-as.matrix(varib.ord[1:50,])
colnames(imp.genes)<-c("MeanDecreaseGini")
write.xlsx(imp.genes,"6.50 Most Important Genes.xlsx",col.names=T,row.names=T)

#Plotting the "imp.genes" matrix
#For visual reasons the order was reversed so as to have an ascending order in the plot
imp.genes.inc<-as.matrix(imp.genes[order(imp.genes,decreasing=F),])
imp.genes.inc.names<-rownames(imp.genes.inc)
plot(y=seq(1:length(imp.genes.inc)),x=imp.genes.inc,pch=19,yaxt="n",xlab="Importance",ylab="", main="50 Most Important Genes",cex.main=1.7,cex.lab=1.3,cex.axis=1.2,cex=1.3,col="steelblue4")
axis(2,at=1:length(imp.genes.inc),labels=imp.genes.inc.names,las=1,cex.axis=0.7)