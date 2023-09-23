library(clusterProfiler)
library(org.Hs.eg.db)
msigdb = clusterProfiler::read.gmt("h.all.v7.4.symbols.gmt")

load("DEG_cancer.Rdata")

#markers<-lapply(markers,function(x){
#	if(length(grep("^RPL|^RPS",rownames(x)))>0){
#		x[-grep("^RPL|^RPS",rownames(x)),]
#	}else{
#		x
#	}
#})
kkk<-sapply(markers,function(x){dim(x)[1]})[c("cMo","ncMo","Macro","DC","pDC")]
pdf("bar_cancer.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]],markers[[5]])
length(unique(markers2$geneName))
hallmark<-list()
hallmarkSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	hallmark[[i]]<-enricher(geneList1, TERM2GENE = msigdb)
	hallmarkSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(hallmark[[i]])$ID,
						p.adj=summary(hallmark[[i]])$p.adjust,
						GeneRatio=eval(parse(text=summary(hallmark[[i]])$GeneRatio)))
}
hallmarkSummary<-rbind(hallmarkSummary_p[[1]],hallmarkSummary_p[[2]],hallmarkSummary_p[[3]],hallmarkSummary_p[[4]],hallmarkSummary_p[[5]])
hallmarkSummary$cellType<-factor(hallmarkSummary$cellType,levels=c("cMo","ncMo","Macro","DC","pDC"))
hallmarkSummary$term<-factor(hallmarkSummary$term,levels=names(sort(table(hallmarkSummary$term))))
library(ggplot2)
pdf("hallmark_cancer.pdf",7,5)
ggplot(hallmarkSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()
dev.off()
write.csv(hallmarkSummary,file="hallmarkSummary_cancer.csv",row.names=F)



load("DEG_normal.Rdata")

#markers<-lapply(markers,function(x){
#	if(length(grep("^RPL|^RPS",rownames(x)))>0){
#		x[-grep("^RPL|^RPS",rownames(x)),]
#	}else{
#		x
#	}
#})
kkk<-sapply(markers,function(x){dim(x)[1]})[c("cMo","ncMo","Macro","DC","pDC")]
pdf("bar_normal.pdf",7,4)
barplot(kkk,col="#8B4513",space=0.8)
dev.off()

markers<-lapply(markers,function(x){x$geneName<-rownames(x);x})
markers2<-rbind(markers[[1]],markers[[2]],markers[[3]],markers[[4]],markers[[5]])
length(unique(markers2$geneName))
hallmark<-list()
hallmarkSummary_p<-list()
for(i in 1:length(markers)){
	geneList1<-rownames(markers[[i]])
	hallmark[[i]]<-enricher(geneList1, TERM2GENE = msigdb)
	hallmarkSummary_p[[i]]<-data.frame(cellType=names(markers)[i],
						term=summary(hallmark[[i]])$ID,
						p.adj=summary(hallmark[[i]])$p.adjust,
						GeneRatio=eval(parse(text=summary(hallmark[[i]])$GeneRatio)))
}
hallmarkSummary<-rbind(hallmarkSummary_p[[1]],hallmarkSummary_p[[2]],hallmarkSummary_p[[3]],hallmarkSummary_p[[4]],hallmarkSummary_p[[5]])
hallmarkSummary$cellType<-factor(hallmarkSummary$cellType,levels=c("cMo","ncMo","Macro","DC","pDC"))
hallmarkSummary$term<-factor(hallmarkSummary$term,levels=names(sort(table(hallmarkSummary$term))))
library(ggplot2)
pdf("hallmark_normal.pdf",7,3)
ggplot(hallmarkSummary,aes(x=cellType,y=term,size=GeneRatio,color=p.adj))+
	geom_point(aes(color=p.adj))+
	scale_color_gradient(low = "#DC143C",high = "black")+
	theme_bw()
dev.off()
write.csv(hallmarkSummary,file="hallmarkSummary_normal.csv",row.names=F)






hallmark2<-lapply(hallmark,summary)
hallmark2<-lapply(hallmark2,function(x){x["HALLMARK_APOPTOSIS","geneID"]})
hallmark2<-lapply(hallmark2,function(x){strsplit(x,"/")[[1]]})
intersect(hallmark2[[1]],intersect(hallmark2[[2]],intersect(hallmark2[[3]],intersect(hallmark2[[4]],hallmark2[[5]]))))

hallmark2<-lapply(hallmark,summary)
hallmark2<-lapply(hallmark2,function(x){x["HALLMARK_HYPOXIA","geneID"]})
hallmark2<-lapply(hallmark2,function(x){strsplit(x,"/")[[1]]})
intersect(hallmark2[[1]],intersect(hallmark2[[2]],intersect(hallmark2[[3]],intersect(hallmark2[[4]],hallmark2[[5]]))))

hallmark2<-lapply(hallmark,summary)
hallmark2<-lapply(hallmark2,function(x){x["HALLMARK_INFLAMMATORY_RESPONSE","geneID"]})
hallmark2<-lapply(hallmark2,function(x){strsplit(x,"/")[[1]]})
intersect(hallmark2[[1]],intersect(hallmark2[[2]],intersect(hallmark2[[3]],intersect(hallmark2[[4]],hallmark2[[5]]))))




hallmark2<-lapply(hallmark,summary)
hallmark2<-lapply(hallmark2,function(x){x["HALLMARK_MYC_TARGETS_V1","geneID"]})
hallmark2<-lapply(hallmark2,function(x){strsplit(x,"/")[[1]]})
intersect(hallmark2[[1]],intersect(hallmark2[[2]],intersect(hallmark2[[3]],intersect(hallmark2[[4]],hallmark2[[5]]))))




