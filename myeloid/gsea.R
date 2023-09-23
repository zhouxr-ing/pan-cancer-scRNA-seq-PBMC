library(clusterProfiler)
library(ggpubr)
options(stringsAsFactors=F)


########################################cancer
library(clusterProfiler)
library(org.Hs.eg.db)
#msigdb = clusterProfiler::read.gmt("../ËèÏµÏ¸°û/h.all.v7.4.symbols.gmt")
msigdb = clusterProfiler::read.gmt("../ËèÏµÏ¸°û/c5.go.bp.v7.5.1.symbols.gmt")

load("DEG.Rdata")
lapply(markers,head)

markers<-lapply(markers,function(x){
	if(length(grep("^RPL|^RPS|^MT-",rownames(x)))>0){
		#x<-x[-grep("^RPL|^RPS|^MT-",rownames(x)),]
	}
	x[which(abs(x$avg_log2FC)>0.25),]
})


for(i in 1:length(markers)){
kkk<-bitr(rownames(markers[[i]]),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db", drop = TRUE)
markers[[i]]$gene<-kkk$ENTREZID[match(rownames(markers[[i]]),kkk$SYMBOL)]
geneList<-markers[[i]]$avg_log2FC
names(geneList)<-rownames(markers[[i]])
geneList<-geneList[order(geneList,decreasing=T)]

library(ggplot2)
library(clusterProfiler)

#geneset <- read.gmt("../ËèÏµÏ¸°û/h.all.v7.4.symbols.gmt")
geneset <- clusterProfiler::read.gmt("../ËèÏµÏ¸°û/c5.go.bp.v7.5.1.symbols.gmt")

#geneset<-geneset[which(geneset[,1]=="HALLMARK_IL6_JAK_STAT3_SIGNALING"),]
egmt <- GSEA(geneList,pvalueCutoff=1,minGSSize=0,TERM2GENE=geneset, verbose=T)
head(egmt)
write.csv(egmt,file=paste0("d:/",names(markers)[i],"_egmt.csv"))
}

p1<-gseaplot(egmt,'GOBP_NEGATIVE_REGULATION_OF_GROWTH',color="#0033CC",color.line="#990033",title = "GLYCOLYSIS") 
p2<-gseaplot(egmt,'GOBP_SPERM_DNA_CONDENSATION',color="#0033CC",color.line="#990033",title = "GOBP_SPERM_DNA_CONDENSATION")
#p3<-gseaplot(egmt,'HALLMARK_ADIPOGENESIS',color="#0033CC",color.line="#990033",title = "ADIPOGENESIS")
#p4<-gseaplot(egmt,'HALLMARK_OXIDATIVE_PHOSPHORYLATION',color="#0033CC",color.line="#990033",title = "OXIDATIVE_PHOSPHORYLATION")
#p5<-gseaplot(egmt,'HALLMARK_ANGIOGENESIS',color="#0033CC",color.line="#990033",title = "ANGIOGENESIS")
p6<-gseaplot(egmt,'HALLMARK_IL2_STAT5_SIGNALING',color="#0033CC",color.line="#990033",title = "IL2_STAT5_SIGNALING")
#p7<-gseaplot(egmt,'HALLMARK_APOPTOSIS',color="#0033CC",color.line="#990033",title = "APOPTOSIS")
p8<-gseaplot(egmt,'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',color="#0033CC",color.line="#990033",title = "EPITHELIAL_MESENCHYMAL_TRANSITION")
p9<-gseaplot(egmt,'HALLMARK_OXIDATIVE_PHOSPHORYLATION',color="#0033CC",color.line="#990033",title = "OXIDATIVE_PHOSPHORYLATION")

pdf("gsea.pdf",14,16)
ggarrange(p1,p2,ncol = 2,nrow =1,widths = c(1,1),heights = c(1,1))
dev.off()

pdf("OXIDATIVE_PHOSPHORYLATION.pdf",7,8)
p9
dev.off()





load("cluster.markers2.Rdata")
DEG<-as.data.frame(cluster.markers)
head(DEG)

DEG$gene<-rownames(DEG)
kkk<-bitr(DEG$gene,"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db", drop = TRUE)
DEG$gene<-kkk$ENTREZID[match(DEG$gene,kkk$SYMBOL)]

geneList<-log2(exp(DEG$avg_logFC))
names(geneList)<-DEG$gene
geneList<-geneList[order(geneList,decreasing=T)]

library(ggplot2)
library(clusterProfiler)

geneset <- read.gmt("h.all.v7.4.entrez.gmt")
geneset<-geneset[which(geneset[,1] %in% 
	paste0("HALLMARK_",c("GLYCOLYSIS","HYPOXIA","ADIPOGENESIS",
				   "OXIDATIVE_PHOSPHORYLATION","ANGIOGENESIS",
				   "IL2_STAT5_SIGNALING","APOPTOSIS",
				   "EPITHELIAL_MESENCHYMAL_TRANSITION",
				   "OXIDATIVE_PHOSPHORYLATION"))),]
egmt <- GSEA(geneList,pvalueCutoff=1,minGSSize=0,TERM2GENE=geneset, verbose=T)
head(egmt)
p1<-gseaplot(egmt,'HALLMARK_GLYCOLYSIS',color="#0033CC",color.line="#990033",title = "GLYCOLYSIS") 
p2<-gseaplot(egmt,'HALLMARK_HYPOXIA',color="#0033CC",color.line="#990033",title = "HYPOXIA")
#p3<-gseaplot(egmt,'HALLMARK_ADIPOGENESIS',color="#0033CC",color.line="#990033",title = "ADIPOGENESIS")
#p4<-gseaplot(egmt,'HALLMARK_OXIDATIVE_PHOSPHORYLATION',color="#0033CC",color.line="#990033",title = "OXIDATIVE_PHOSPHORYLATION")
#p5<-gseaplot(egmt,'HALLMARK_ANGIOGENESIS',color="#0033CC",color.line="#990033",title = "ANGIOGENESIS")
p6<-gseaplot(egmt,'HALLMARK_IL2_STAT5_SIGNALING',color="#0033CC",color.line="#990033",title = "IL2_STAT5_SIGNALING")
#p7<-gseaplot(egmt,'HALLMARK_APOPTOSIS',color="#0033CC",color.line="#990033",title = "APOPTOSIS")
p8<-gseaplot(egmt,'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',color="#0033CC",color.line="#990033",title = "EPITHELIAL_MESENCHYMAL_TRANSITION")
p9<-gseaplot(egmt,'HALLMARK_OXIDATIVE_PHOSPHORYLATION',color="#0033CC",color.line="#990033",title = "OXIDATIVE_PHOSPHORYLATION")

pdf("gsea2.pdf",14,16)
ggarrange(p1,p2,p6,p8,ncol = 2,nrow =2,widths = c(1,1),heights = c(1,1))
dev.off()

pdf("OXIDATIVE_PHOSPHORYLATION2.pdf",7,8)
p9
dev.off()






DEG<-read.delim("Layer02_logNor_FCcluster2.0_DEG.txt",stringsAsFactors=F)
DEG<-DEG[which(DEG$cluster==2),]
DEG<-DEG[which(DEG$p_val_adj<0.05),]
dim(DEG)
head(DEG)

kkk<-bitr(DEG$gene,"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db", drop = TRUE)
DEG$gene<-kkk$ENTREZID[match(DEG$gene,kkk$SYMBOL)]
#geneList<-DEG$avg_logFC
geneList<-DEG$pct.1
names(geneList)<-DEG$gene
geneList<-geneList[order(geneList,decreasing=T)]

enrich <- enricher(names(geneList),TERM2GENE=geneset); head(enrich)
dotplot(enrich)











