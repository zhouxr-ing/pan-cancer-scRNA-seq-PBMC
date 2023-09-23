library(clusterProfiler)
library(org.Hs.eg.db)

load("markers_identify.Rdata")
markers<-markers_identify[which(markers_identify$avg_log2FC<1),]
geneList1<-rownames(markers)
geneList1<-geneList1[-grep("^RPL|^RPS",geneList1)]
eg1 = bitr(geneList1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ekk1 <- enrichKEGG(gene=eg1$ENTREZID, keyType="kegg", organism='hsa', qvalueCutoff = 0.05) #KEGG富集分析
print(head(ekk1))
print(barplot(ekk1))


geneList1 <- sort(geneList1, decreasing=T)
geneList2 <- sort(geneList2, decreasing=T)

eg1 = bitr(names(geneList1), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ekk1 <- enrichKEGG(gene=eg1$ENTREZID, keyType="kegg", organism='hsa', qvalueCutoff = 0.05) #KEGG富集分析
egoBP1 <- enrichGO(gene = geneList1, keyType="SYMBOL", OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff = 0.05, readable= F) #GO_BP
print(head(egoBP1))
print(dotplot(egoBP1))
egoMF1 <- enrichGO(gene = geneList1, keyType="SYMBOL", OrgDb="org.Hs.eg.db", ont="MF", pvalueCutoff = 0.05, readable= F) #GO_MF
print(head(egoMF1))
dotplot(egoMF1)

