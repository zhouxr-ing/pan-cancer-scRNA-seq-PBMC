########################################cancer
library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../髓系细胞/h.all.v7.4.symbols.gmt")

load("DEG.Rdata")
markers<-markers[which(names(markers)!="low_quality")]
lapply(markers,head)

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC>0.25 & x$p_val<0.05),]
})

kkk<-sapply(markers,function(x){dim(x)[1]})[c("CD4_Tn","CD4_Tm","Treg",
                                              "CD8_Tn","GZMB_CD8","GZMK_CD8","cell_cycle","MAIT","low_quality")]

pdf("bar_cancer.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]],
		markers[[5]],markers[[6]],markers[[7]],markers[[8]])
length(unique(markers2$geneName))

BP<-list()
BPSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	BP[[i]]<-enrichGO(geneList1,'org.Hs.eg.db',keyType='SYMBOL',ont='BP')
	BPSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(BP[[i]])$Description,
						p.adj=summary(BP[[i]])$p.adjust,
						GeneRatio=sapply(summary(BP[[i]])$GeneRatio,function(x){eval(parse(text=x))}),
						geneID=summary(BP[[i]])$geneID
)
}
BPSummary<-rbind(BPSummary_p[[1]],BPSummary_p[[2]])
if(length(BPSummary_p)>2){
  for(i in 3:length(BPSummary_p)){
    BPSummary<-rbind(BPSummary,BPSummary_p[[i]])
  }
}
BPSummary$cellType<-factor(BPSummary$cellType,levels=c("CD4_Tn","CD4_Tm","Treg",
                                                         "CD8_Tn","GZMB_CD8","GZMK_CD8","cell_cycle","MAIT"))
BPSummary$term<-factor(BPSummary$term,levels=names(sort(table(BPSummary$term))))
write.csv(BPSummary,file="BPSummary_cancer.csv",row.names=F)
intersectedTerm<-c("negative regulation of hemopoiesis",
		"positive regulation of leukocyte differentiation",
		"positive regulation of mRNA metabolic process",
		"negative regulation of lymphocyte apoptotic process",
		"positive regulation of T cell activation",
		"positive regulation of response to endoplasmic reticulum stress",
		"positive regulation of intracellular transport",
		"negative regulation of leukocyte apoptotic process",
		"positive regulation of pri-miRNA transcription by RNA polymerase II",
		"positive regulation of intracellular protein transport",
		"negative regulation of T cell apoptotic process",
		"positive regulation of nuclear-transcribed mRNA catabolic process, deadenylation-dependent decay",
		"positive regulation of leukocyte activation",
		"positive regulation of RNA splicing",
		"positive regulation of I-kappaB kinase/NF-kappaB signaling",
		"positive regulation of cell activation",
		"positive regulation of cellular catabolic process",
		"positive regulation of T cell differentiation",
		"negative regulation of interleukin-2 production",
		"negative regulation of receptor signaling pathway via JAK-STAT",
		"negative regulation of viral transcription",
		"negative regulation of viral genome replication",
		"negative regulation of leukocyte activation",
		"negative regulation of interleukin-2 production",
		"positive regulation of interleukin-6 biosynthetic process",
		"positive regulation of chemokine production")

library(ggplot2)
pdf("BP_cancer.pdf",9.5,5)
BPSummary<-BPSummary[which(BPSummary$term %in% intersectedTerm),]
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

BPSummary1<-BPSummary


########################################normal

library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../髓系细胞/h.all.v7.4.symbols.gmt")

load("DEG.Rdata")
markers<-markers[which(names(markers)!="low_quality")]

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC<(-0.25) & x$p_val<0.05),]
})

kkk<-sapply(markers,function(x){dim(x)[1]})[c("CD4_Tn","CD4_Tm","Treg",
                                              "CD8_Tn","GZMB_CD8","GZMK_CD8",
							"cell_cycle","MAIT")]
pdf("bar_normal.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]],
		markers[[5]],markers[[6]],markers[[7]],markers[[8]])
length(unique(markers2$geneName))

BP<-list()
BPSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	BP[[i]]<-enrichGO(geneList1,'org.Hs.eg.db',keyType='SYMBOL',ont='BP',pAdjustMethod="none",pvalueCutoff=1,qvalueCutoff=1)
	if(dim(summary(BP[[i]]))[1]>0){
	  BPSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(BP[[i]])$Description,
						pvalue=summary(BP[[i]])$pvalue,
						GeneRatio=sapply(summary(BP[[i]])$GeneRatio,function(x){eval(parse(text=x))}),
						geneID=summary(BP[[i]])$geneID
						)
	}
}
BPSummary<-rbind(BPSummary_p[[1]],BPSummary_p[[2]])
if(length(BPSummary_p)>2){
  for(i in 3:length(BPSummary_p)){
    BPSummary<-rbind(BPSummary,BPSummary_p[[i]])
  }
}

BPSummary$cellType<-factor(BPSummary$cellType,levels=c("CD4_Tn","CD4_Tm","Treg",
                                                       "CD8_Tn","GZMB_CD8","GZMK_CD8","cell_cycle","MAIT","low_quality"))
BPSummary$term<-factor(BPSummary$term,levels=names(sort(table(BPSummary$term))))
BPSummary<-BPSummary[which(BPSummary$pvalue<0.05),]

intersectedTerm<-c("positive regulation of protein oligomerization",
			"regulation of protein oligomerization",
			"negative regulation of cytokine production",
			"negative regulation of exocytosis",
			"response to tumor necrosis factor",
			"DNA deamination",
			"tumor necrosis factor-mediated signaling pathway",
			"double-strand break repair",
			"regulation of DNA metabolic process",
			"DNA integrity checkpoint",
			"cell cycle DNA replication",
			"negative regulation of DNA metabolic process")
library(ggplot2)
pdf("BP_normal.pdf",10,4)
#BPSummary<-BPSummary[which(BPSummary$term %in% intersectedTerm),]
#BPSummary$term<-factor(BPSummary$term,levels=intersectedTerm)
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

BPSummary1$cellType<-paste0(BPSummary1$cellType,"_c")
BPSummary2$cellType<-paste0(BPSummary2$cellType,"_n")
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
