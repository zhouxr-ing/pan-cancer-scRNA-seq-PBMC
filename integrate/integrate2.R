load("result/mergeGEO/train.Rdata")

library(Seurat)
cat("\n\n----------train set Seurat pipeline----------\n")
pbmc.list<-list()
for(i in 1:length(train)){
  pbmc.list[[i]]<-CreateSeuratObject(train[[i]],min.cells=3)
  pbmc.list[[i]]<-SCTransform(pbmc.list[[i]],verbose=F)
  print(i)
}

pbmcMerge<-merge(pbmc.list[[1]],y=c(pbmc.list[[2]],pbmc.list[[3]],pbmc.list[[4]],pbmc.list[[5]],pbmc.list[[6]],pbmc.list[[7]],pbmc.list[[8]]))
pbmcMerge<-FindVariableFeatures(pbmcMerge,nfeatures=dim(pbmcMerge)[1]) 
pbmcMerge<-ScaleData(pbmcMerge)
length(VariableFeatures(pbmcMerge)) 
pbmcMerge<-RunPCA(pbmcMerge)
pbmcMerge<-RunTSNE(pbmcMerge,check_duplicates=F)

pdf("QC.pdf",height=5,width=12)
VlnPlot(pbmcMerge,features=c("nFeature_RNA","nCount_RNA"),ncol=2,pt.size=0)
dev.off()

pdf("result/mergeGEO/before_train.pdf")
DimPlot(object=pbmcMerge,reduction="tsne")
dev.off()
print("000")

pbmc.features<-SelectIntegrationFeatures(object.list=pbmc.list,nfeatures=3000)
print("111")
pbmc.list<-PrepSCTIntegration(object.list=pbmc.list,anchor.features=pbmc.features)
print("222")
pbmc.anchors<-FindIntegrationAnchors(object.list=pbmc.list,dims=1:20,k.anchor=5,k.filter=30,normalization.method="SCT",anchor.features=pbmc.features)
print("333")
pbmc.integrated<-IntegrateData(anchorset = pbmc.anchors,normalization.method = "SCT")
print("444")
save(pbmc.integrated,file="result/mergeGEO/pbmc.integrated_train.Rdata")

DefaultAssay(pbmc.integrated)<-"integrated"
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
pbmc.integrated<-RunPCA(pbmc.integrated)
pbmc.integrated<-RunTSNE(pbmc.integrated,check_duplicates=F)

pdf("result/mergeGEO/after_train.pdf")
DimPlot(pbmc.integrated,reduction="tsne",group.by='orig.ident')
dev.off()






rm(list=ls())
load("result/mergeGEO/test.Rdata")
cat("\n\n----------test set Seurat pipeline----------\n")
pbmc.list<-list()
for(i in 1:length(test)){
  pbmc.list[[i]]<-CreateSeuratObject(test[[i]],min.cells=3)
  pbmc.list[[i]]<-SCTransform(pbmc.list[[i]],verbose = FALSE)
  print(i)
}
pbmcMerge<-merge(pbmc.list[[1]],y=c(pbmc.list[[2]],pbmc.list[[3]],pbmc.list[[4]],pbmc.list[[5]],pbmc.list[[6]],pbmc.list[[7]],pbmc.list[[8]]))
pbmcMerge<-FindVariableFeatures(pbmcMerge,nfeatures=dim(pbmcMerge)[1])
pbmcMerge<-ScaleData(pbmcMerge)
length(VariableFeatures(pbmcMerge))
pbmcMerge<-RunPCA(pbmcMerge)
pbmcMerge<-RunTSNE(pbmcMerge,check_duplicates=F)

pdf("QC.pdf",height=5,width=12)
VlnPlot(pbmcMerge,features=c("nFeature_RNA","nCount_RNA"),ncol=2,pt.size=0)
dev.off()

pdf("result/mergeGEO/before_test.pdf")
DimPlot(object=pbmcMerge,reduction="tsne")
dev.off()

print("000")
pbmc.features<-SelectIntegrationFeatures(object.list=pbmc.list,nfeatures=3000)
print("111")
pbmc.list<-PrepSCTIntegration(object.list=pbmc.list,anchor.features=pbmc.features)
print("222")
pbmc.anchors<-FindIntegrationAnchors(object.list=pbmc.list,dims=1:20,k.anchor=5,k.filter=30,normalization.method="SCT",anchor.features=pbmc.features)
print("333")
pbmc.integrated<-IntegrateData(anchorset = pbmc.anchors,normalization.method = "SCT")
print("444")
save(pbmc.integrated,file="result/mergeGEO/pbmc.integrated.test.Rdata")

DefaultAssay(pbmc.integrated)<-"integrated"
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
pbmc.integrated<-RunPCA(pbmc.integrated)
pbmc.integrated<-RunTSNE(pbmc.integrated,check_duplicates=F)

pdf("result/mergeGEO/after_test.pdf")
DimPlot(pbmc.integrated,reduction="tsne",group.by='orig.ident')
dev.off()




