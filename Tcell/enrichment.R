library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../髓系细胞/h.all.v7.4.symbols.gmt")

load("DEG.Rdata")

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC>0.25 & x$p_val_adj<0.05),]
})
names(markers)<-c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle","Treg")
kkk<-sapply(markers,function(x){dim(x)[1]})[c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle","Treg")]
pdf("bar_cancer.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],
			markers[[4]],markers[[5]],markers[[6]],
			markers[[7]],markers[[8]],markers[[9]])

length(unique(markers2$geneName))
hallmark<-list()
hallmarkSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	hallmark[[i]]<-enricher(geneList1, TERM2GENE = msigdb)
	hallmarkSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(hallmark[[i]])$ID,
						p.adj=summary(hallmark[[i]])$p.adjust,
						GeneRatio=sapply(summary(hallmark[[i]])$GeneRatio,function(x){eval(parse(text=x))})
)
}
hallmarkSummary<-rbind(hallmarkSummary_p[[1]],hallmarkSummary_p[[2]],hallmarkSummary_p[[3]],
				hallmarkSummary_p[[4]],hallmarkSummary_p[[5]],hallmarkSummary_p[[6]],
				hallmarkSummary_p[[7]],hallmarkSummary_p[[8]],hallmarkSummary_p[[9]])
hallmarkSummary$cellType<-factor(hallmarkSummary$cellType,levels=c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle","Treg"))
hallmarkSummary$term<-factor(hallmarkSummary$term,levels=names(sort(table(hallmarkSummary$term))))
library(ggplot2)
pdf("hallmark_cancer.pdf",7,5)
ggplot(hallmarkSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
write.csv(hallmarkSummary,file="hallmarkSummary_cancer.csv",row.names=F)

########################################normal

library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../髓系细胞/h.all.v7.4.symbols.gmt")

load("DEG.Rdata")


markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC<(-0.25) & x$p_val_adj<0.05),]
})
names(markers)<-c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle","Treg")
kkk<-sapply(markers,function(x){dim(x)[1]})[c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle","Treg")]
pdf("bar_normal.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

kkk<-sapply(markers,function(x){dim(x)[1]})[c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle","Treg")]
pdf("bar_cancer.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],
			markers[[4]],markers[[5]],markers[[6]],
			markers[[7]],markers[[8]],markers[[9]])

length(unique(markers2$geneName))
hallmark<-list()
hallmarkSummary_p<-list()
for(i in 1:8){############Treg富集不上
	geneList1<-rownames(markers[[i]])
	hallmark[[i]]<-enricher(geneList1, TERM2GENE = msigdb)
	hallmarkSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(hallmark[[i]])$ID,
						p.adj=summary(hallmark[[i]])$p.adjust,
						GeneRatio=eval(parse(text=summary(hallmark[[i]])$GeneRatio)))
}
hallmarkSummary<-rbind(hallmarkSummary_p[[1]],hallmarkSummary_p[[2]],hallmarkSummary_p[[3]],
				hallmarkSummary_p[[4]],hallmarkSummary_p[[5]],hallmarkSummary_p[[6]],
				hallmarkSummary_p[[7]],hallmarkSummary_p[[8]])############Treg富集不上，去掉hallmarkSummary_p[[9]]
hallmarkSummary$cellType<-factor(hallmarkSummary$cellType,levels=c("CD4","GZMB_CD8","naive_CD4","KLRB1_CD8","naive_CD8","GZMK_CD8","low_quality","cell_cycle"))
hallmarkSummary$term<-factor(hallmarkSummary$term,levels=names(sort(table(hallmarkSummary$term))))
library(ggplot2)
pdf("hallmark_normal.pdf",7,5)
ggplot(hallmarkSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()+
	theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
write.csv(hallmarkSummary,file="hallmarkSummary_normal.csv",row.names=F)


