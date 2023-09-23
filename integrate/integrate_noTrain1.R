allData<-list()
allData$GSE117988<-read.delim("../pretreatment/matrix/GSE117988_PBMC_naive.txt",stringsAsFactors=F,row.names=1)
allData$GSE123192<-read.delim("../pretreatment/matrix/GSE123192.txt",stringsAsFactors=F,row.names=1)
allData$GSE145281muc<-read.delim("../pretreatment/matrix/GSE145281_PBMC_mUC.txt",stringsAsFactors=F,row.names=1)
allData$GSE145281rcc<-read.delim("../pretreatment/matrix/GSE145281_PBMC_RCC.txt",stringsAsFactors=F,row.names=1)
allData$GSE148190<-read.delim("../pretreatment/matrix/GSE148190_PBMC.txt",stringsAsFactors=F,row.names=1)
allData$GSE139324<-read.delim("../pretreatment/matrix/GSE139324_PBMC.txt",stringsAsFactors=F,row.names=1)
allData$normal1<-read.delim("../pretreatment/matrix/GSE139324_normal2.txt",stringsAsFactors=F,row.names=1)
allData$normal2<-read.delim("../pretreatment/matrix/GSE128066_PBMC.txt",stringsAsFactors=F,row.names=1)

cat("----------allData----------\n")
for(i in 1:length(allData)){
  print(dim(allData[[i]]))
  print(allData[[i]][1:4,1:4])
}

cat("\n\n----------min.features=200----------\n")
allData<-lapply(allData,function(a){a[,which(apply(a,2,function(x){length(which(x>0))})>200)]})
for(i in 1:length(allData)){
  print(dim(allData[[i]]))
}

#allDataing set
set.seed(123)
allData$GSE117988<-allData$GSE117988[,sample(1:dim(allData$GSE117988)[2],round(dim(allData$GSE117988)[2]/5))]
allData$GSE123192<-allData$GSE123192[,sample(1:dim(allData$GSE123192)[2],round(dim(allData$GSE123192)[2]/10))]
allData$GSE145281muc<-allData$GSE145281muc[,sample(1:dim(allData$GSE145281muc)[2],round(dim(allData$GSE145281muc)[2]/5))]
allData$GSE145281rcc<-allData$GSE145281rcc[,sample(1:dim(allData$GSE145281rcc)[2],round(dim(allData$GSE145281rcc)[2]/5))]
allData$GSE148190<-allData$GSE148190[,sample(1:dim(allData$GSE148190)[2],round(dim(allData$GSE148190)[2]/5))]
allData$GSE139324<-allData$GSE139324[,sample(1:dim(allData$GSE139324)[2],round(dim(allData$GSE139324)[2]/5))]
allData$normal1<-allData$normal1[,sample(1:dim(allData$normal1)[2],round(dim(allData$normal1)[2]/2.5))]
allData$normal2<-allData$normal2[,sample(1:dim(allData$normal2)[2],round(dim(allData$normal2)[2]/2.5))]

#"MCC","BC","mUC","RCC","MM","normal1"
names(allData$GSE117988)<-paste("MCC",names(allData$GSE117988),sep="_")
names(allData$GSE123192)<-paste("BC",names(allData$GSE123192),sep="_")
names(allData$GSE145281muc)<-paste("UC",names(allData$GSE145281muc),sep="_")
names(allData$GSE145281rcc)<-paste("RCC",names(allData$GSE145281rcc),sep="_")
names(allData$GSE148190)<-paste("MM",names(allData$GSE148190),sep="_")
names(allData$GSE139324)<-paste("HNSCC",names(allData$GSE139324),sep="_")
names(allData$normal1)<-paste("normal1",names(allData$normal1),sep="_")
names(allData$normal2)<-paste("normal2",names(allData$normal2),sep="_")
save(allData,file="result/mergeGEO/allData.Rdata")

cat("\n\n----------allData----------\n")
for(i in 1:length(allData)){
  print(dim(allData[[i]]))
}

library(Seurat)
cat("\n\n----------allData set Seurat pipeline----------\n")
pbmc.list<-list()
for(i in 1:length(allData)){
  pbmc.list[[i]]<-CreateSeuratObject(allData[[i]])
  pbmc.list[[i]]<-SCTransform(pbmc.list[[i]],verbose = FALSE)
  print(i)
}
print("000")

pbmc.features<-SelectIntegrationFeatures(object.list=pbmc.list,nfeatures=3000)
print("111")
pbmc.list<-PrepSCTIntegration(object.list=pbmc.list,anchor.features=pbmc.features)
print("222")
pbmc.anchors<-FindIntegrationAnchors(object.list=pbmc.list,dims=1:20,k.anchor=5,k.filter=30,normalization.method="SCT",anchor.features=pbmc.features)
print("333")
pbmc.integrated<-IntegrateData(anchorset = pbmc.anchors,normalization.method = "SCT")
print("444")
save(pbmc.integrated,file="result/mergeGEO/pbmc.integrated_allData.Rdata")

pbmc.integrated<-RunPCA(pbmc.integrated)
pbmc.integrated<-RunTSNE(pbmc.integrated,check_duplicates=F)

pdf("result/mergeGEO/integrated_tsne_allData.pdf")
DimPlot(pbmc.integrated,reduction="tsne",group.by='orig.ident')
dev.off()


