library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../ËèÏµÏ¸°û/h.all.v7.4.symbols.gmt")

load("DEG.Rdata")

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC>0.25),]
})
length(unique(c(rownames(markers[[1]]),rownames(markers[[2]]),rownames(markers[[3]]))))

names(markers)<-c("naive_memory_B","plasma_B","immaturate_B")
kkk<-sapply(markers,function(x){dim(x)[1]})[c("immaturate_B","naive_memory_B","plasma_B")]
pdf("bar_cancer.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]])

length(unique(markers2$geneName))
hallmark<-list()
hallmarkSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	hallmark[[i]]<-enricher(geneList1, TERM2GENE = msigdb)
	hallmarkSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(hallmark[[i]])$ID,
						p.adj=summary(hallmark[[i]])$p.adjust,
						GeneRatio=sapply(summary(hallmark[[i]])$GeneRatio,function(x){eval(parse(text=x))}),
						geneID=summary(hallmark[[i]])$geneID)
}
hallmarkSummary<-rbind(hallmarkSummary_p[[1]],hallmarkSummary_p[[2]],hallmarkSummary_p[[3]])
hallmarkSummary$cellType<-factor(hallmarkSummary$cellType,levels=c("immaturate_B","naive_memory_B","plasma_B"))
hallmarkSummary$term<-factor(hallmarkSummary$term,levels=names(sort(table(hallmarkSummary$term))))
library(ggplot2)
pdf("hallmark_cancer.pdf",7,5)
ggplot(hallmarkSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()
dev.off()
write.csv(hallmarkSummary,file="hallmarkSummary_cancer.csv",row.names=F)

########################################normal

library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("../ËèÏµÏ¸°û/h.all.v7.4.symbols.gmt")

load("DEG.Rdata")
markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS",rownames(x)))>0){
		x<-x[-grep("^RPL|^RPS",rownames(x)),]
	}
	x[which(x$avg_log2FC<(-0.25)),]
})
length(unique(c(rownames(markers[[1]]),rownames(markers[[2]]),rownames(markers[[3]]))))

names(markers)<-c("naive_memory_B","plasma_B","immaturate_B")
kkk<-sapply(markers,function(x){dim(x)[1]})[c("immaturate_B","naive_memory_B","plasma_B")]
pdf("bar_normal.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]])

length(unique(markers2$geneName))
hallmark<-list()
hallmarkSummary_p<-list()
for(i in 2:length(markers)){
	geneList1<-rownames(markers[[i]])
	hallmark[[i]]<-enricher(geneList1, TERM2GENE = msigdb,minGSSize=0)
	hallmarkSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(hallmark[[i]])$ID,
						p.adj=summary(hallmark[[i]])$p.adjust,
						GeneRatio=eval(parse(text=summary(hallmark[[i]])$GeneRatio)),
						geneID=summary(hallmark[[i]])$geneID)
}
hallmarkSummary<-rbind(hallmarkSummary_p[[1]],hallmarkSummary_p[[2]],hallmarkSummary_p[[3]])
hallmarkSummary$cellType<-factor(hallmarkSummary$cellType,levels=c("immaturate_B","naive_memory_B","plasma_B"))
hallmarkSummary$term<-factor(hallmarkSummary$term,levels=names(sort(table(hallmarkSummary$term))))
library(ggplot2)
pdf("hallmark_normal.pdf",7,5)
ggplot(hallmarkSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()
dev.off()
write.csv(hallmarkSummary,file="hallmarkSummary_normal.csv",row.names=F)

kkk<-lapply(hallmarkSummary[which(hallmarkSummary[,2]=="HALLMARK_MYC_TARGETS_V1"),"geneID"],function(x){strsplit(x,"/")[[1]]})
intersect(kkk[[1]],kkk[[2]])


