library(pheatmap)
options(stringsAsFactors=FALSE)
load("result/cluster/pbmc.integrated.cluster.Rdata")
DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
geneMat<-pbmc.integrated@assays$RNA@scale.data

load("RRA_featureName.Rdata")
featureName

pdf("RRA.pdf",height=5,width=30)
pheatmap(geneMat[which(rownames(geneMat) %in% featureName),])
dev.off()













library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

options(stringsAsFactors=FALSE)
load("result/cluster/pbmc.integrated.cluster.Rdata")
DefaultAssay(pbmc.integrated)<-"RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
pbmc.integrated<-FindVariableFeatures(pbmc.integrated,nfeatures=dim(pbmc.integrated)[1])
pbmc.integrated<-ScaleData(pbmc.integrated)
geneMat<-pbmc.integrated@assays$RNA@data

load("RRA_featureName.Rdata")
geneAnn<-featureName

cellAnn<-pbmc.integrated@meta.data$orig.ident
cellAnn[grep("normal",cellAnn)]<-"normal"
cellAnn[which(cellAnn!="normal")]<-"tumor"
cellAnn<-data.frame(sample_type=cellAnn)

colors = function(n, type){
  nmax = brewer.pal.info[type, 1]
  colors = c()
  if(n < 3){
    colors = brewer.pal(3, type)[1:n]
  }else if(n <= nmax){
    colors = brewer.pal(n, type)
  }else{
    colors = colorRampPalette(brewer.pal(nmax, type), space="Lab")(n)
  }
  colors
}

geneMat<-geneMat[which(rownames(geneMat) %in% geneAnn),]

mycol <- colorRampPalette(c("white", "firebrick3"))(100)
col1 <- list(sample_type=c(tumor="pink", normal="forestgreen"))

ha0<-HeatmapAnnotation(df=cellAnn,col=col1)
mat<-as.matrix(geneMat)
mat<-t(mat)
ht0<-Heatmap(mat, name = "log2(TPM+1)", col = mycol,
        split=factor(cellAnn$sample_type),
   #     column_title = "Cluster mehtod: ward.D, distance: euclidean", 
        row_names_gp = gpar(fontsize = 4), column_names_gp = gpar(fontsize = 4),
        show_row_names = T,
row_names_side="left",
     #   cluster_columns = T,
        clustering_method_rows = "ward.D",
        column_title_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(color_bar = "continuous"))


pdf("exprHeatmap_kk.pdf",12,12)
ht_list = ht0
draw(ht_list)
dev.off()
