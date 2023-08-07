# ROC

rm(list = ls()) 
options(stringsAsFactors = F)

if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library(pacman)
p_load(data.table,pROC,limma)

biof<-read.csv("gse7084_2507.csv")
biof=as.matrix(biof)
rownames(biof)=biof[,1]
exp=biof[,-1]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

sample=read.table("clinical7084_2507.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)


group <- c(rep("0",8),rep("1",7)) 
group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])


biof2<-merge(y, x,by = "row.names", all = T)
colnames(biof2)[1]="id"
colnames(biof2)[2]="status"
Veendata<-read.table("hubgene7084.txt",header=F,sep="\t",check.names = F)
veenExp=biof2[,c("id","status",as.vector(Veendata[,1]))]


for(wenyun in colnames(veenExp[,3:ncol(veenExp)])){
  myroc=roc(veenExp$status, veenExp[,wenyun])
  myci=ci.auc(myroc, method="bootstrap")
  mycinum=as.numeric(myci)
  pdf(file=paste0(wenyun,"_","ROC",".pdf"))   
  plot(myroc, print.auc=T, col="#00AFBB", legacy.axes=T, main=wenyun)   
  text(0.40, 0.45, paste0("95% CI: ",sprintf("%.03f",mycinum[1]),"-",sprintf("%.03f",mycinum[3])), col="#00AFBB")
  dev.off()
}

