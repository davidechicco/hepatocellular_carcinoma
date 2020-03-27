# Enviroment
setwd(".")
rm(list=ls())
set.seed(13)
cat("\014")

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
  ris = c(ris,caret::negPredValue(data=as.factor(YF),reference=as.factor(YT),negative="0"))
}

# Results
method = "DT"
D <- read.csv('Data.csv',header=TRUE,sep = ',')
if (method == "DT") {
  v = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,28,29,50) }
if (method == "SVL" || method == "SVK" || method == "MLP") {
  v = c(50) }
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
  if (method == "DT") {
    strmethod = "rpart2"
    grid = expand.grid(maxdepth=c(2,4,6,8,10,12,14)) }
  if (method == "SVL") {
    strmethod = "svmLinear"
    grid = expand.grid(C=c(.0001,.0005,.001,.005,.01,.05,.1,.5,1,5,10,50)) }
  if (method == "SVK") {
    strmethod = "svmRadial"
    grid = expand.grid(C=c(.0001,.0005,.001,.005,.01,.05,.1,.5,1,5,10,50),
                       sigma=c(.0001,.0005,.001,.005,.01,.05,.1,.5,1,5,10,50)) }
  if (method == "MLP") {
    strmethod = "mlpKerasDropout"; 
    grid = expand.grid(size=c(5,10,20,40,80,160),
                       dropout=c(0,.001,.01,.1),           
                       batch_size=c(nl/10,nl),
                       lr=c(.001,.01,.1,1),
                       rho=c(.9,0.09),
                       decay=c(.001,.01,.1,1),activation=c("relu")) }
  trctrl = trainControl(method="repeatedcv",
                        number=10,
                        repeats=1,
                        sampling="up",
                        allowParallel=TRUE)
  M <- train(x=XL,
             y=YL,
             method=strmethod,
             trControl=trctrl,
             tuneGrid=grid,
             preProcess=c("center","scale"),
             metric="Kappa")
  YF = predict(M,XT)
  YYT = c(YYT,as.numeric(YT)-1)
  YYF = c(YYF,as.numeric(YF)-1)
}
cat("\014")
#
ris = caret::confusionMatrix(as.factor(YYF), as.factor(YYT),positive="1")
print(ris[2])
ris = c()
for (i in c(1:mc)) {
  j = sample(length(YYF),replace=TRUE)
  ris = rbind(ris, MetricsClass(YYF[j],YYT[j])) }
s = sprintf("mcc:    %.3f pm %.3f",mean(ris[,1]),sd(ris[,1])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=FALSE,sep="\n")
s = sprintf("f1:     %.3f pm %.3f",mean(ris[,2]),sd(ris[,2])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("acc:    %.3f pm %.3f",mean(ris[,3]),sd(ris[,3])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("rec:    %.3f pm %.3f",mean(ris[,4]),sd(ris[,4])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("spec:   %.3f pm %.3f",mean(ris[,5]),sd(ris[,5])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("prauc:  %.3f pm %.3f",mean(ris[,6]),sd(ris[,6])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("rocauc: %.3f pm %.3f",mean(ris[,7]),sd(ris[,7])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("ppv:    %.3f pm %.3f",mean(ris[,8]),sd(ris[,8])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")
s = sprintf("ppm:    %.3f pm %.3f",mean(ris[,9]),sd(ris[,9])); print(s); cat(s,file=sprintf("RIS/%s_Met.txt",method),append=TRUE,sep="\n")