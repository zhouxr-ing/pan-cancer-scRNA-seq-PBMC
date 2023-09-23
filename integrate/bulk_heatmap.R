options(stringsAsFactors=FALSE)

library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

library(limma)
geneMat<-read.delim("data/bulk/GSE120691_ccvsn.txt")
geneMat2<-as.matrix(geneMat[,-(1:2)])
rownames(geneMat2)<-geneMat[,2]
geneMat<-avereps(geneMat2)
geneMat2<-apply(geneMat,2,function(x){log2(as.numeric(x)+1)})
rownames(geneMat2)<-rownames(geneMat)
geneMat<-geneMat2

load("RRA_featureName.Rdata")
featureName
geneAnn<-data.frame(regulation=c(rep("up_regulation",11),rep("down_regulation",20)))
rownames(geneAnn)<-featureName

sampleAnn<-colnames(geneMat)
sampleAnn[grep("NC",sampleAnn)]<-"normal"
sampleAnn[which(sampleAnn!="normal")]<-"tumor"
sampleAnn<-data.frame(sample_type=sampleAnn)
rownames(sampleName)<-colnames(geneMat)

geneMat<-geneMat[which(rownames(geneMat) %in% rownames(geneAnn)),]

mycol <- colorRampPalette(c("white", "firebrick3"))(100)
col0 <- list(regulation=c(up_regulation="pink", down_regulation="forestgreen"))
col1 <- list(sample_type=c(tumor="#FFD700", normal="#87CEFA"))

ha0<-rowAnnotation(df=sampleAnn,annotation_legend_param = list(ncol = 1, title = "Type", title_position = "topcenter"),
                          width = unit(4, "mm"), col = col1)
ha_gene<-HeatmapAnnotation(geneAnn,col = col0)

mat<-normalizeBetweenArrays(t(as.matrix(geneMat)))
ht0<-Heatmap(mat, name = "log2(TPM+1)", col = mycol,
        split=factor(sampleAnn$sample_type),
   #     column_title = "Cluster mehtod: ward.D, distance: euclidean", 
        row_names_gp = gpar(fontsize = 4), column_names_gp = gpar(fontsize = 4),
        show_row_names = T,
row_names_side="left",
        cluster_columns = F,
        clustering_method_rows = "ward.D",
        column_title_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(color_bar = "continuous"),
        top_annotation = ha_gene)


pdf("result/bulk/exprHeatmap_kk.pdf",12,12)
ht_list = ht0+ha0
draw(ht_list)
dev.off()

