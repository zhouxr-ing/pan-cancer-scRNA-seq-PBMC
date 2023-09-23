library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("h.all.v7.4.symbols.gmt")

load("DEG_cancer.Rdata")

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC>0.25 & x$p_val<0.05),]
})
kkk<-sapply(markers,function(x){dim(x)[1]})[c("cMo","ncMo","Macro","DC","pDC")]
pdf("bar_cancer.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]],markers[[5]])

length(unique(markers2$geneName))
BP<-list()
BPSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	BP[[i]]<-enrichGO(geneList1,'org.Hs.eg.db',keyType='SYMBOL',ont='BP')
	BPSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(BP[[i]])$Description,
						pvalue=summary(BP[[i]])$pvalue,
						GeneRatio=sapply(summary(BP[[i]])$GeneRatio,function(x){eval(parse(text=x))}),
						geneID=summary(BP[[i]])$geneID
)
}
BP1<-BP
BPSummary<-rbind(BPSummary_p[[1]],BPSummary_p[[2]],BPSummary_p[[3]],BPSummary_p[[4]],BPSummary_p[[5]])
BPSummary$cellType<-factor(BPSummary$cellType,levels=c("cMo","ncMo","Macro","DC","pDC"))
BPSummary$term<-factor(BPSummary$term,levels=names(sort(table(BPSummary$term))))

library(ggplot2)
pdf("BP_cancer.pdf",7,5)
BPSummary<-BPSummary[which(BPSummary$pvalue<0.05),]
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=pvalue))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()
dev.off()
write.csv(BPSummary,file="BPSummary_cancer.csv",row.names=F)
BPSummary1<-BPSummary

intersectedTerm<-c("positive regulation of cytokine production",
"positive regulation of cell activation",
"negative regulation of intrinsic apoptotic signaling pathway",
"positive regulation of lymphocyte proliferation",
"negative regulation of apoptotic signaling pathway",
"positive regulation of leukocyte activation",
"positive regulation of angiogenesis",
"positive regulation of response to biotic stimulus",
"positive regulation of interleukin-6 production",
"negative regulation of intrinsic apoptotic signaling pathway",
"positive regulation of apoptotic signaling pathway",
"positive regulation of neuron apoptotic process",
"positive regulation of tumor necrosis factor production",
"positive regulation of epithelial cell proliferation",
"positive regulation of interferon-beta production",
"positive regulation of myeloid cell differentiation",
"positive regulation of type I interferon production",
"negative regulation of cytokine secretion",
"positive regulation of cell activation",
"negative regulation of apoptotic signaling pathway",
"positive regulation of leukocyte activation",
"positive regulation of angiogenesis",
"positive regulation of response to biotic stimulus",
"positive regulation of interleukin-6 production",
"positive regulation of inflammatory response",
"negative regulation of intrinsic apoptotic signaling pathway",
"positive regulation of apoptotic signaling pathway",
"positive regulation of neuron apoptotic process",
"positive regulation of tumor necrosis factor production",
"positive regulation of epithelial cell proliferation",
"positive regulation of interferon-beta production",
"negative regulation of immune system process",
"negative regulation of protein kinase activity",
"positive regulation of myeloid cell differentiation",
"positive regulation of type I interferon production",
"negative regulation of cell activation",
"negative regulation of leukocyte activation",
"negative regulation of lymphocyte activation",
"negative regulation of cytokine secretion",
"antigen processing and presentation of exogenous peptide antigen via MHC class II",
"antigen processing and presentation of peptide antigen via MHC class II",
"antigen processing and presentation of peptide or polysaccharide antigen via MHC class II")


library(ggplot2)
pdf("BP_cancer.pdf",9.5,5)
BPSummary<-BPSummary[which(BPSummary$term %in% unique(intersectedTerm)),]
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=pvalue))+
	geom_point(aes(color=pvalue))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

########################################normal

library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("h.all.v7.4.symbols.gmt")

load("DEG_normal.Rdata")

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC>0.25 & x$p_val<0.05),]
})

kkk<-sapply(markers,function(x){dim(x)[1]})[c("cMo","ncMo","Macro","DC","pDC")]
pdf("bar_normal.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]],markers[[5]])

length(unique(markers2$geneName))
BP<-list()
BPSummary_p<-list()
for(i in c(1:5)){
	geneList1<-rownames(markers[[i]])
	BP[[i]]<-enrichGO(geneList1,'org.Hs.eg.db',keyType='SYMBOL',ont='BP',pvalueCutoff=1,qvalueCutoff=1)
	BPSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(BP[[i]])$Description,
						pvalue=summary(BP[[i]])$pvalue,
						GeneRatio=sapply(summary(BP[[i]])$GeneRatio,function(x){eval(parse(text=x))}),
						geneID=summary(BP[[i]])$geneID
)
}
BP2<-BP
BPSummary<-rbind(BPSummary_p[[1]],BPSummary_p[[2]],BPSummary_p[[3]],BPSummary_p[[4]],BPSummary_p[[5]])
BPSummary$cellType<-factor(BPSummary$cellType,levels=c("cMo","ncMo","Macro","DC","pDC"))
BPSummary$term<-factor(BPSummary$term,levels=names(sort(table(BPSummary$term))))
BPSummary<-BPSummary[which(BPSummary$pvalue<0.05),]


intersectedTerm2<-c("positive regulation of cell activation",
"negative regulation of intrinsic apoptotic signaling pathway",
"positive regulation of lymphocyte proliferation",
"negative regulation of apoptotic signaling pathway",
"positive regulation of leukocyte activation",
"positive regulation of angiogenesis",
"positive regulation of response to biotic stimulus",
"positive regulation of interleukin-6 production",
"negative regulation of intrinsic apoptotic signaling pathway",
"positive regulation of apoptotic signaling pathway",
"positive regulation of neuron apoptotic process",
"positive regulation of tumor necrosis factor production",
"positive regulation of epithelial cell proliferation",
"positive regulation of interferon-beta production",
"positive regulation of myeloid cell differentiation",
"positive regulation of type I interferon production",
"negative regulation of cytokine secretion",
"positive regulation of cell activation",
"negative regulation of apoptotic signaling pathway",
"positive regulation of leukocyte activation",
"positive regulation of angiogenesis",
"positive regulation of response to biotic stimulus",
"positive regulation of interleukin-6 production",
"positive regulation of inflammatory response",
"negative regulation of intrinsic apoptotic signaling pathway",
"positive regulation of apoptotic signaling pathway",
"positive regulation of neuron apoptotic process",
"positive regulation of tumor necrosis factor production",
"positive regulation of epithelial cell proliferation",
"positive regulation of interferon-beta production",
"negative regulation of immune system process",
"negative regulation of protein kinase activity",
"positive regulation of myeloid cell differentiation",
"positive regulation of type I interferon production",
"negative regulation of cell activation",
"negative regulation of leukocyte activation",
"negative regulation of lymphocyte activation",
"negative regulation of cytokine secretion",
"antigen processing and presentation of exogenous peptide antigen via MHC class II",
"antigen processing and presentation of peptide antigen via MHC class II",
"antigen processing and presentation of peptide or polysaccharide antigen via MHC class II")

library(ggplot2)
pdf("BP_normal.pdf",9,4)
BPSummary<-BPSummary[which(BPSummary$term %in% unique(intersectedTerm2)),]
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=pvalue))+
	geom_point(aes(color=pvalue))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()
dev.off()
write.csv(BPSummary,file="BPSummary_normal.csv",row.names=F)
BPSummary2<-BPSummary

ppp<-union(strsplit(summary(BP[[3]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]],
		strsplit(summary(BP[[2]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]])
ppp3<-strsplit(summary(BP[[3]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]]
ppp2<-strsplit(summary(BP[[2]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]]

kkk<-unique(c(strsplit(summary(BP1[[3]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]],
		strsplit(summary(BP1[[2]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]],
		strsplit(summary(BP1[[1]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]]))
kkk3<-strsplit(summary(BP1[[3]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]]
kkk2<-strsplit(summary(BP1[[2]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]]
kkk1<-strsplit(summary(BP1[[3]])["HALLMARK_MYC_TARGETS_V1","geneID"][1],"/")[[1]]

BPSummary1$cellType<-paste0(as.character(BPSummary1$cellType),"_c")
BPSummary2$cellType<-paste0(as.character(BPSummary2$cellType),"_n")
BPSummary<-rbind(BPSummary1,BPSummary2)
BPSummary<-BPSummary[which(BPSummary$GeneRatio>0.05),]
ppp<-reshape2::acast(BPSummary,term~cellType)
ppp[which(is.na(ppp))]<-0
kkk<-pheatmap::pheatmap(ppp)

BPSummary$term<-factor(BPSummary$term,levels=kkk$tree_row$labels[kkk$tree_row$order])
pdf("BP.pdf",15,14)
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

png("BP.png",1000,1000)
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
