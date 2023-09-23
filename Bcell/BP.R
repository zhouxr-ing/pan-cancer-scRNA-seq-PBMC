library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../髓系细胞/h.all.v7.4.symbols.gmt")

load("E:/硕士毕设/GEOdata_jd/GEOdata/Integrate4/3.cellType/Bcell/DEG.Rdata")

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC>0.25),]
})
names(markers)<-c("naive_B","memory_B","plasma_B","immaturate_B")
kkk1<-sapply(markers,function(x){dim(x)[1]})[c("immaturate_B","naive_B","memory_B","plasma_B")]
names(kkk1)<-paste0(names(kkk1),"_c")
pdf("bar_cancer.pdf",7,4)
barplot(kkk1,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]])

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
BPSummary<-rbind(BPSummary_p[[1]],BPSummary_p[[2]],BPSummary_p[[3]],BPSummary_p[[4]])
BPSummary$cellType<-factor(BPSummary$cellType,levels=c("immaturate_B","naive_B","memory_B","plasma_B"))
BPSummary$term<-factor(BPSummary$term,levels=names(sort(table(BPSummary$term))))

library(ggplot2)
pdf("BP_cancer.pdf",10,5)
BPSummary<-BPSummary[which(BPSummary$p.adj<0.05),]
#PSummary<-BPSummary[which(BPSummary$term %in% termList),]
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
write.csv(BPSummary,file="BPSummary_cancer.csv",row.names=F)
BPSummary1<-BPSummary


########################################normal

library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../髓系细胞/h.all.v7.4.symbols.gmt")

load("E:/硕士毕设/GEOdata_jd/GEOdata/Integrate4/3.cellType/Bcell/DEG.Rdata")

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC<(-0.25)),]
})
names(markers)<-c("naive_B","memory_B","plasma_B","immaturate_B")
kkk2<-sapply(markers,function(x){dim(x)[1]})[c("immaturate_B","naive_B","memory_B","plasma_B")]
names(kkk2)<-paste0(names(kkk2),"_n")
kkk<-kkk2

pdf("bar_normal.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]])

length(unique(markers2$geneName))
BP<-list()
BPSummary_p<-list()
for(i in 3:length(markers)){
	geneList1<-rownames(markers[[i]])
	BP[[i]]<-enrichGO(geneList1,'org.Hs.eg.db',keyType='SYMBOL',ont='BP',pAdjustMethod="none")
	BPSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(BP[[i]])$Description,
						p.adj=summary(BP[[i]])$p.adjust,
						GeneRatio=sapply(summary(BP[[i]])$GeneRatio,function(x){eval(parse(text=x))}),
						geneID=summary(BP[[i]])$geneID
)
}
BPSummary<-rbind(BPSummary_p[[1]],BPSummary_p[[2]],BPSummary_p[[3]],BPSummary_p[[4]])
BPSummary$cellType<-factor(BPSummary$cellType,levels=c("immaturate_B","naive_B","memory_B","plasma_B"))
BPSummary$term<-factor(BPSummary$term,levels=names(sort(table(BPSummary$term))))
BPSummary<-BPSummary[which(BPSummary$p.adj<0.05),]

library(ggplot2)
pdf("BP_normal.pdf",12,5)
#PSummary<-BPSummary[which(BPSummary$term %in% termList),]
ggplot(BPSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
write.csv(BPSummary,file="BPSummary_normal.csv",row.names=F)
BPSummary2<-BPSummary


