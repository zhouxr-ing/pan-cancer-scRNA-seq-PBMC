library(Seurat)

load("result/cluster/pbmc.integrated.cluster.Rdata")
DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)

pbmc.integrated@meta.data$label<-pbmc.integrated@meta.data$orig.ident
pbmc.integrated@meta.data$label[which(pbmc.integrated@meta.data$label %in% c("normal1","normal2"))]<-"normal"
pbmc.integrated@meta.data$label[which(pbmc.integrated@meta.data$label!="normal")]<-"cancer"

marker0<-c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4")#inhibitor receptor
pbmc.box<-GetAssayData(object=pbmc.integrated,slot='data')
pbmc.box<-pbmc.box[which(rownames(pbmc.box) %in% marker0),]
pbmc.box<-data.frame(t(pbmc.box),label=pbmc.integrated@meta.data$label)


library(reshape2)
pbmc.box<-melt(pbmc.box,data.id.vars=c("label"),variable.name="gene")
pbmc.box[,3]<-as.numeric(pbmc.box[,3])

library(ggplot2)
pdf("kkk.pdf")
ggplot(data=pbmc.box,aes(gene,value,fill=label))+
	theme_bw()+
	geom_violin()+
	scale_fill_brewer()
dev.off()

