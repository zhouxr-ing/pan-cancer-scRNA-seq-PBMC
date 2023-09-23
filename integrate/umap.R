library(Seurat)
load("result/mergeGEO/pbmc.integrated_train.Rdata")

DefaultAssay(pbmc.integrated)<-"integrated"
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
pbmc.integrated<-RunPCA(pbmc.integrated)
pbmc.integrated<-RunUMAP(pbmc.integrated,dims=1:15)

pdf("result/mergeGEO/after_train_umap.pdf")
DimPlot(pbmc.integrated,reduction="umap",group.by='orig.ident')
dev.off()






library(Seurat)
cat("----------train set----------\n")
setwd("/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/GEOdata/Integrate/")
load("result/mergeGEO/pbmc.integrated_train.Rdata")

pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
pbmc.integrated<-RunPCA(pbmc.integrated)
pbmc.integrated<-RunTSNE(pbmc.integrated,check_duplicates=F)

cat("\n\n----------find clusters----------\n")
pbmc.integrated<-FindNeighbors(pbmc.integrated,dims = 1:15)
pbmc.integrated<-FindClusters(pbmc.integrated,resolution=0.5)
pbmc.integrated<-RunTSNE(pbmc.integrated,dims=1:15,check_duplicates=F)
pbmc.integrated<-RunUMAP(pbmc.integrated,dims=1:15)



pdf("result/cluster/pbmc_integrated_cluster.pdf",15,4)
plot1<-DimPlot(pbmc.integrated,reduction="tsne",group.by="orig.ident")
plot2<-DimPlot(pbmc.integrated,reduction="tsne",label=T)
plot3<-DimPlot(pbmc.integrated,reduction="umap",label=T)
CombinePlots(plots=list(plot1,plot2,plot3),ncol=3)
pbmc.integrated<-RunUMAP(pbmc.integrated,dims=1:15)

dev.off()
#save(pbmc.integrated,file="result/cluster/pbmc.integrated.cluster.Rdata")



DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)

rownames(pbmc.integrated@assays$RNA@data)[grep("^CD16",rownames(pbmc.integrated@assays$RNA@data))]
which(rownames(pbmc.integrated@assays$RNA@data)=="CD4")

marker0<-c("TCF7","SELL","LEF1","CCR7",#naive
	"LAG3","TIGIT","PDCD1","HAVCR2","CTLA4",#inhibitor receptor
	"IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7",#cytokiness
	"MS4A1","CD79A",#B
        "CD3D","CD3E","CD3G",#T
	"CD8A","CD8B",
        "FOXP3","IL2RA",#Treg
        "NKG7","CST3",#NK
        "CD14","CST3","LYZ",#macrophage, monocyte, neutrophil
	"LYZ","FCER1A"#DC
)

pdf("result/cluster/marker0.pdf",20,20)
p<-DoHeatmap(pbmc.integrated,features=marker0)
print(p)
dev.off()

pdf("result/cluster/marker_kk1.pdf",20,24)
p<-FeaturePlot(object=pbmc.integrated,features=marker0[1:16])
print(p)
dev.off()

pdf("result/cluster/marker_kk2.pdf",20,24)
p<-FeaturePlot(object=pbmc.integrated,features=marker0[17:32])
print(p)
dev.off()



Heatmap(scep,features=top10$gene,group.by="label")

marker0<-list(cluster7=c("MS4A7","FCGR3A","ITIFM3","SPN","B2M"),
        cluster11=c("HMGB2","TUBA1B","HSP90B1","PPIB","IGHA1","JCHAIN","STMN1","ITM2C","MZB1","IGHA2","GAPDH"),
        naive=c("TCF7","SELL","LEF1","CCR7","CD8A","CD8B"),#naive
        inhibitor=c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4"),#inhibitor receptor
        cytokiness=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),#cytokiness
        Bcell=c("MS4A1","CD79A"),#B
        Tcell=c("CD3D","CD3E","CD3G"),#T
        CD8=c("CD8A","CD8B","CD4"),
        Treg=c("FOXP3","IL2RA"),#Treg
        NK=c("NKG7","CST3"),#NK
        myeloid=c("CD14","CST3","LYZ","FCGR3B"),#macrophage, monocyte, neutrophil
        DC=c("LYZ","FCER1A")#DC
)

for(i in 1:length(marker0)){
	width_pdf<-(sign(floor(length(marker0[[i]])/4))-1)*2.5*(4-length(marker0[[i]]))+10
	height_pdf<-ceiling(length(marker0[[i]])/4)*2
	pdf(paste0("result/cluster/marker_",names(marker0)[i],"_kk.pdf"),width_pdf,height_pdf)
	p<-FeaturePlot(object=pbmc.integrated,features=marker0[[i]],reduction="tsne",ncol=width_pdf/2.5)
	print(p)
	dev.off()

}



"HMGB2","TUBA1B","HSP90B1","PPIB","IGHA1","JCHAIN","STMN1","ITM2C","MZB1","IGHA2"


 
