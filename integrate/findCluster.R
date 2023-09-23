library(Seurat)
cat("----------train set----------\n")
load("result/mergeGEO/pbmc.integrated_train.Rdata")

pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
pbmc.integrated<-RunPCA(pbmc.integrated)
pbmc.integrated<-RunTSNE(pbmc.integrated,check_duplicates=F)

cat("\n\n----------find clusters----------\n")
pbmc.integrated<-FindNeighbors(pbmc.integrated,dims = 1:15)
pbmc.integrated<-FindClusters(pbmc.integrated,resolution=0.5)
pbmc.integrated<-RunTSNE(pbmc.integrated,dims=1:15,check_duplicates=F)

pdf("result/cluster/pbmc_integrated_cluster.pdf",15,4)
plot1<-DimPlot(pbmc.integrated,reduction="tsne",group.by="orig.ident")
plot2<-DimPlot(pbmc.integrated,reduction="tsne",label=T)
plot3<-DimPlot(pbmc.integrated,reduction="tsne",group.by="orig.ident",cells.highlight=rownames(pbmc.integrated@meta.data)[which(pbmc.integrated@meta.data$orig])
CombinePlots(plots=list(plot1,plot2,plot3),ncol=3)
dev.off()
#save(pbmc.integrated,file="result/cluster/pbmc.integrated.cluster.Rdata")

data<-as.data.frame(t(as.matrix(pbmc.integrated@assays$RNA@counts)))
data$cell_type<-pbmc.integrated@meta.data$seurat_cluster
data$group<-pbmc.integrated@meta.data$orig.ident
data$group[which(data$group %in% c("normal1","normal2"))]<-"normal"
data$group[which(data$group!="normal")]<-"cancer"
save(data,file="iTALK_data.Rdata")

pdf("result/cluster/normal_tumor.pdf")
DimPlot(pbmc.integrated,reduction="tsne",label=T,
        cells.highlight=rownames(pbmc.integrated@meta.data[which(pbmc.integrated@meta.data$orig.ident %in% c("normal1","normal2")),]),
        cols.highlight="lightblue")
dev.off()

cat("\n\n----------marker heatmap----------\n")
library(magrittr)
library(dplyr)
DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
pbmc.markers<-FindAllMarkers(pbmc.integrated,only.pos=T,
                          logfc.threshold=0.5,
		     min.diff.pct=0.3,
                     min.pct=0)
save(pbmc.markers,file="result/cluster/pbmc.markers.Rdata")
top10<-pbmc.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)
pdf("result/cluster/markerHeatmap2.pdf",20,15)
p<-DoHeatmap(pbmc.integrated,features=top10$gene)
print(p)
dev.off()

cat("\n\n----------marker plot----------\n")
table(pbmc.integrated@meta.data$seurat_clusters)
for(i in as.character(unique(pbmc.markers[,6]))){
#for(i in c(1:8,10:12)){
        eachTop<-top10[which(top10[,6]==i),7]
        png(paste0("result/cluster/feature2_",i,".png"),1000,750)
        p<-FeaturePlot(object=pbmc.integrated,features=c(eachTop[[1]]))
        print(p)
        dev.off()
#        pdf(paste0("result/cluster/feature_",i,".pdf"),10,7.5)
#        print(p)
#        dev.off()
	print(i)
}









marker0<-c('HBA1','CD3G','CD4','CD8A','CD8B','IL7R','CCR7','SELL','CCR6','KLRB1','CTLA4','TIGIT','LAG3','FOXP3','IL2RA','NR4A1','TRAC','TRDC','TRGC1','FCGR3A','CD160','GNLY','NKG7','KLRF1','GZMH','GZMA','PRF1','NKG7','XCL1','BCL11A','TCL1A','CD123','MS4A1','JCHAIN','MZB1','SSR4','CD14','LYZ','CD1C','CD1E','FCER1A','CLEC9A','CLEC10A','XCR1','PPBP','CPA3','MKI67','CD14','ITGA2B','VGFR','NRP1')

marker0<-c("MNDA","CSTA","GPX1","FCN1","CD14","LYZ","VCAN","S100A12","LYZ","S100A8","S100A9","TRAC","CD3D","IL32","GIMAP7","LTB","PIK3IP1","MAL","TCF7","IL7R","LEF1","NOSIP","LDHB","CCR7","HLA−DPB1","TCL1A","LINC00926","IGKC","CD79B","IGHM","HLA−DRA","HLA−DQA1","MS4A1","CD79A","KLRF1","SPON2","GZMK","GZMH","COTL1","SAT1","SERPINA1","CFD","FCER1G","MS4A7","AIF1","FCGR3A","IFITM3","LST1","GABARAPL1","SELENOK","RNF125","TNFAIP3","TMEM2","FAM177A1","LEPROTL1","RNF19A","CTD−2192J16.22","CD34","CD41A","CD61")


pdf("kkk.pdf")
#pdf("result/cluster/markerHeatmap_kk.pdf",15,20)
p<-DoHeatmap(pbmc.integrated,features=marker0)
print(p)
dev.off()

png(paste0("result/cluster/featurekk_",1,".png"),1500,2000)
p<-FeaturePlot(object=pbmc.integrated,features=marker0)
print(p)
dev.off()


##QC
pdf("result/QC_cluster.pdf",height=5,width=12)
VlnPlot(pbmc.integrated,features=c("nFeature_RNA","nCount_RNA"),ncol=2,pt.size=0,group.by="seurat_clusters")
dev.off()




###
load("result/cluster/pbmc.integrated.cluster.Rdata")
DefaultAssay(pbmc.integrated)<-"RNA"

pbmc.integrated@meta.data$myCluster<-pbmc.integrated@meta.data$seurat_clusters
clusterLabel<-data.frame(seurat_cluster=0:12,cell_type=c("monocyte","naive_Tcell","naive_Tcell","Bcell","NKcell","CD8_Tcell_GZMKhigh","CD8_Tcell_GZMBhigh","mature myeloid cell","un-known","Tcell9","DCcell","cycle cell","cycle plasma cell"))
pbmc.integrated@meta.data$myCluster<-clusterLabel[match(pbmc.integrated@meta.data$seurat_clusters,clusterLabel[,1]),2]
pbmc.integrated@meta.data$label<-pbmc.integrated@meta.data$orig.ident
pbmc.integrated@meta.data$label[which(pbmc.integrated@meta.data$label %in% c("normal1","normal2"))]<-"normal"
pbmc.integrated@meta.data$label[which(pbmc.integrated@meta.data$label!="normal")]<-"cancer"

#bar plot
pdf("cancer_normal_bar.pdf",8,3)
ggplot(data = pbmc.integrated@meta.data, aes(x = myCluster, fill=label))+geom_bar(stat = 'count',position = 'fill')+labs(x = '')+coord_flip() 
dev.off()



pdf("result/cluster/pbmc_cluster_label.pdf",7,4)
DimPlot(pbmc.integrated,reduction="tsne",group.by="myCluster")
dev.off()

library(magrittr)
library(dplyr)
clusterLabel2<-unique(clusterLabel[,2])
for(i in 1:length(clusterLabel2)){
    scep<-subset(pbmc.integrated,cells=rownames(pbmc.integrated@meta.data[which(pbmc.integrated@meta.data$myCluster==clusterLabel2[i]),]))
    
    scep<-FindNeighbors(scep,dims = 1:15)
    scep<-FindClusters(scep,resolution=0.5)
    scep<-RunTSNE(scep,dims=1:15,check_duplicates=F,perplexity=30)

    pdf(paste0("result/cluster/cluster_p30_part",clusterLabel2[i],".pdf"),15,4)
    plot1<-DimPlot(scep,reduction="tsne",group.by="label")
    plot2<-DimPlot(scep,reduction="tsne",group.by="orig.ident")
    plot3<-DimPlot(scep,reduction="tsne",label=T)
    p<-CombinePlots(plots=list(plot1,plot2,plot3),ncol=3)
    print(p)
    dev.off()

    DefaultAssay(scep)<-"RNA"
    scep<-NormalizeData(scep, verbose = FALSE)
    scep<-FindVariableFeatures(scep,nfeatures=dim(scep)[1])
    scep<-ScaleData(scep)
    pbmc.markers<-FindAllMarkers(scep,only.pos=T)

    top10<-pbmc.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)
    pdf(paste0("result/cluster/markerHeat_part",i,".pdf"),20,20)
    p<-DoHeatmap(scep,features=top10$gene)
    print(p)
    dev.off()
}



DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)

pdf("result/cluster/kkk.pdf")
DoHeatmap(scep,features=top10$gene,group.by="label")
dev.off()









######cluster0_7之间的marker
marker0_7<-FindMarkers(object = pbmc.integrated, ident.1 = 0, ident.2 = 7)
marker0_7<-c("MS4A7","FCGR3A","ITIFM3","SPN","B2M")

DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)

pdf("result/cluster/heatmap0_7.pdf",6,6)
DoHeatmap(pbmc.integrated,features=rownames(marker0_7)[1:20],group.by="seurat_clusters")
dev.off()

marker3_11<-FindMarkers(object = pbmc.integrated, ident.1 = 3, ident.2 = 11)
head(x = markers)

DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)

pdf("result/cluster/marker0_7.pdf",10,4)
p<-FeaturePlot(object=pbmc.integrated,features=marker0_7)
print(p)
dev.off()


pdf("result/cluster/heatmap3_11.pdf",6,6)

DoHeatmap(subset(pbmc.integrated,cells=rownames(pbmc.integrated@meta.data[which(pbmc.integrated@meta.data$seurat_clusters %in% c(3,11)),])),
	features=c(rownames(marker3_11)[1:20],"MS4A1","CD79A"),group.by="seurat_clusters")
dev.off()


marker3_11_con<-FindConservedMarkers(pbmc.integrated, ident.1=3, ident.2=11,grouping.var="seurat_clusters")
marker3_11_con<-FindConservedMarkers(subset(pbmc.integrated,cells=rownames(pbmc.integrated@meta.data[which(pbmc.integrated@meta.data$seurat_clusters %in% c(3,11)),])), ident.1=3, ident.2=11,grouping.var="seurat_clusters")
pdf("result/cluster/heatmap3_11_con.pdf",4,6)
DoHeatmap(subset(subset(pbmc.integrated,cells=rownames(pbmc.integrated@meta.data[which(pbmc.integrated@meta.data$seurat_clusters %in% c(3,11)),])),
        cells=rownames(pbmc.integrated@meta.data[which(pbmc.integrated@meta.data$seurat_clusters %in% c(3,11)),])),
        features=c(rownames(marker3_11_con)[1:20],"MS4A1","CD79A"),group.by="seurat_clusters")
dev.off()




