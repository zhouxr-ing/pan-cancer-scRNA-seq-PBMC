library(limma)
library(Seurat)
logFoldChange=0.5
adjustP=0.05

load("result/cluster/pbmc.integrated.cluster.Rdata")
#load("result/mergeGEO/pbmc.integrated_train.Rdata")
load("result/cluster/pbmc.markers.Rdata")
DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
geneMat<-pbmc.integrated@assays$RNA@scale.data

#marker gene
library(dplyr)
library(magrittr)
top10<-pbmc.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)

feature_svm<-c()
diffGene<-list()
for(i in 0:12){
  geneMat2<-geneMat[,which(pbmc.integrated@meta.data$seurat_clusters==i)]
  class<-pbmc.integrated@meta.data$orig.ident[which(pbmc.integrated@meta.data$seurat_clusters==i)]
  class[grep("normal",class)]<-"normal"
  class[which(class!="normal")]<-"tumor"

  design<-model.matrix(~0+factor(class))
  colnames(design)<-c("con","case")
  fit<-lmFit(geneMat2,design)
  cont.matrix<-makeContrasts(case-con,levels=design)
  fit2<-contrasts.fit(fit,cont.matrix)
  fit2<-eBayes(fit2)

  allDiff<-topTable(fit2,adjust='fdr',number=200000)
  allDiff$adj.P.Val<-allDiff$adj.P.Val+1e-300

  marker10<-top10[which(top10$cluster==i),7][[1]]
  png(paste0("result/DEG/vol_",i,".png"))
  #pdf(paste0("result/DEG/vol_",i,".pdf"))
  yMax<-max(-log10(allDiff$adj.P.Val))
  xMax<-max(abs(allDiff$logFC))

  plot(allDiff$logFC,-log10(allDiff$adj.P.Val),ylab="-log10(adj.P.Val)",xlab="logFC",col="gray",
     main=paste0("cluster",i), ylim=c(0,yMax),xlim=c(-xMax,xMax),xaxs="i",pch=20, cex=0.7)
  diffSub1<-subset(allDiff,adj.P.Val<adjustP & logFC>logFoldChange)
  diffSub2<-subset(allDiff,adj.P.Val<adjustP & logFC<(-logFoldChange))
  diffSub3<-allDiff[which(rownames(allDiff) %in% marker10),]
  diffSub4<-diffSub3[which(rownames(diffSub3) %in% c(rownames(diffSub1),rownames(diffSub2))),]
  diffGene[[i+1]]<-allDiff

  diffSub<-subset(allDiff,adj.P.Val<adjustP & logFC>logFoldChange)
  points(diffSub$logFC,-log10(diffSub$adj.P.Val),pch=16,col=rgb(247,68,97,100,maxColorValue=255),cex=0.9)#red
  diffSub<-subset(allDiff,adj.P.Val<adjustP & logFC<(-logFoldChange))
  points(diffSub$logFC,-log10(diffSub$adj.P.Val),pch=16,col=rgb(160,191,124,100,maxColorValue=255),cex=0.9)#green
  diffSub<-subset(allDiff,adj.P.Val<adjustP & logFC>1)
  points(diffSub$logFC,-log10(diffSub$adj.P.Val),pch=20,col=rgb(255,0,0,100,maxColorValue=255),cex=1.2)#deep red
  diffSub<-subset(allDiff,adj.P.Val<adjustP & logFC<(-1))
  points(diffSub$logFC,-log10(diffSub$adj.P.Val),pch=20,col=rgb(101,147,74,100,maxColorValue=255),cex=1.2)#deep green

  points(diffSub3$logFC,-log10(diffSub3$adj.P.Val),pch=1,col="black",cex=0.7,lwd=1.9)
  
  if(dim(diffSub4)[1]!=0){
    text(x=diffSub4[,1]+0.1,y=-log10(diffSub4[,5]),labels=row.names(diffSub4),cex=0.8)
    feature_svm<-c(feature_svm,row.names(diffSub4))
  }
  abline(h=-log10(adjustP),lty=2,lwd=1.5)
  abline(v=0.5,lty=2,lwd=1.5)
  abline(v=-0.5,lty=2,lwd=1.5)
  abline(v=1,lty=2,lwd=1.5)
  abline(v=-1,lty=2,lwd=1.5)

  dev.off()
  print(i)
}
names(diffGene)<-0:12
save(diffGene,file="result/DEG/diffGene.Rdata")
save(feature_svm,file="result/DEG/feature_svm.Rdata")
