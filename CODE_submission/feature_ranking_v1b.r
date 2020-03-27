# Enviroment
setwd(".")
rm(list=ls())
cat("\014")
set.seed(13)


options(repos = list(CRAN="http://cran.rstudio.com/"))
# update R packages
list.of.packages <- c("easypackages", "randomForest","caret","mltools","MLmetrics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("easypackages")
libraries(list.of.packages)

# Metrics
MetricsClass <- function(YF,YT){
  ris = mltools::mcc(preds=YF,actuals=YT)
  ris = c(ris,MLmetrics::F1_Score(y_true=YT,y_pred=YF,positive=1))
  ris = c(ris,MLmetrics::Accuracy(y_true=YT,y_pred=YF))
  ris = c(ris,MLmetrics::Recall(y_true=YT,y_pred=YF,positive=1))
  ris = c(ris,MLmetrics::Specificity(y_true=YT,y_pred=YF,positive=1))
  ris = c(ris,MLmetrics::PRAUC(y_true=YT,y_pred=YF))
  ris = c(ris,MLmetrics::AUC(y_true=YT,y_pred=YF))
  ris = c(ris,caret::posPredValue(data=as.factor(YF),reference=as.factor(YT),positive="1"))
  ris = c(ris,caret::negPredValue(data=as.factor(YF),reference=as.factor(YT),negative="0")) }

# Results
D <- read.csv('Data.csv',header=TRUE,sep = ',')
v = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,28,29,50)
ix = c(1:49)
iy = c(50)
for (i in v) {
  D[,i] = as.factor(D[,i]) }
#
n = nrow(D)
nl = round(.97*n)
YYT = c()
YYF = c()
mc = 1000
IM1 = array(0,dim=c(length(ix),mc))
IM2 = array(0,dim=c(length(ix),mc))
for (i in c(1:mc)) {
  print(sprintf("mc: %04d",i))
  k = sample(n)
  il = k[1:nl]; it = k[nl:n]
  XL = D[il,ix]; YL = D[il,iy]
  XT = D[it,ix]; YT = D[it,iy]
  tmp = min(sum(YL==levels(YL)[1]),sum(YL==levels(YL)[2]))
  M = randomForest(x=XL,y=YL,
                   mtry = round(sqrt(ncol(XL))),
                   ntree=5000,
                   do.trace=FALSE,
                   importance=TRUE,
                   sampsize=c(tmp,tmp))
  tmp = c(1:length(ix))
  IM1[sort(importance(M,type=1,scale=FALSE),decreasing=TRUE,index.return=TRUE)$ix,i] = tmp 
  IM2[sort(importance(M,type=2,scale=FALSE),decreasing=TRUE,index.return=TRUE)$ix,i] = tmp 
  YF = predict(M,XT)
  YYT = c(YYT,as.numeric(YT)-1)
  YYF = c(YYF,as.numeric(YF)-1) }
cat("\014")
#
I = array(0,dim=length(ix))
for (i in c(1:length(ix))) {
  I[i] = sum(IM1[i,])/mc }
tmp = sort(I,decreasing=FALSE,index.return=TRUE)$ix
cat("",file="RIS/RF_MDA.txt",append=FALSE,sep="")
print(sprintf("MDA"))
for (i in tmp){
  s = sprintf("%02.1f - %s",I[i],names(D)[ix[i]])
  print(s)
  cat(s,file="RIS/RF_MDA.txt",append=TRUE,sep="\n") }
I = array(0,dim=length(ix))
for (i in c(1:length(ix))) {
  I[i] = sum(IM2[i,])/mc }
tmp = sort(I,decreasing=FALSE,index.return=TRUE)$ix
cat("",file="RIS/RF_MDI.txt",append=FALSE,sep="")
print(sprintf("MDI"))
for (i in tmp){
  s = sprintf("%02.1f - %s",I[i],names(D)[ix[i]])
  print(s)
  cat(s,file="RIS/RF_MDI.txt",append=TRUE,sep="\n") }
#
XL = D[,ix]; YL = D[,iy]
tmp = min(sum(YL==levels(YL)[1]),sum(YL==levels(YL)[2]))
FinalModel = randomForest(x=XL,y=YL,
                          mtry = round(sqrt(ncol(XL))),
                          ntree=5000,
                          do.trace=FALSE,
                          importance=TRUE,
                          sampsize=c(tmp,tmp))
par(mfrow=c(1,2))
varImpPlot(FinalModel,type=1,scale=FALSE)
varImpPlot(FinalModel,type=2,scale=FALSE)
write.table(importance(M, type=1, scale=FALSE),
            file=sprintf('RIS/RF_MDA_Round.txt',mc,s), col.names=T, row.names=T)
write.table(importance(M, type=2, scale=FALSE),
            file=sprintf('RIS/RF_MDI_Round.txt',mc,s), col.names=T, row.names=T)
#
ris = caret::confusionMatrix(as.factor(YYF), as.factor(YYT),positive="1")
print(ris[2])
ris = c()
for (i in c(1:mc)) {
  j = sample(length(YYF),replace=TRUE)
  ris = rbind(ris, MetricsClass(YYF[j],YYT[j])) }
s = sprintf("mcc:    %.3f pm %.3f",mean(ris[,1]),sd(ris[,1])); print(s); cat(s,file="RIS/RF_Met.txt",append=FALSE,sep="\n")
s = sprintf("f1:     %.3f pm %.3f",mean(ris[,2]),sd(ris[,2])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("acc:    %.3f pm %.3f",mean(ris[,3]),sd(ris[,3])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("rec:    %.3f pm %.3f",mean(ris[,4]),sd(ris[,4])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("spec:   %.3f pm %.3f",mean(ris[,5]),sd(ris[,5])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("prauc:  %.3f pm %.3f",mean(ris[,6]),sd(ris[,6])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("rocauc: %.3f pm %.3f",mean(ris[,7]),sd(ris[,7])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("ppv:    %.3f pm %.3f",mean(ris[,8]),sd(ris[,8])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")
s = sprintf("ppm:    %.3f pm %.3f",mean(ris[,9]),sd(ris[,9])); print(s); cat(s,file="RIS/RF_Met.txt",append=TRUE,sep="\n")