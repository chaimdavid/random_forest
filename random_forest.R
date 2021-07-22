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
