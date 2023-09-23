data_all<-read.delim("GSE99nnn/GSE99254/suppl/GSE99254_NSCLC.TCell.S12346.count.txt.gz",check.names=F,stringsAsFactors=F)
data_all[1:6,1:6]
data_all<-data_all[-which(is.na(data_all[,2])),]
data_all2<-as.matrix(data_all[,-(1:2)])
rownames(data_all2)<-data_all[,2]
data_pbmc<-data_all2[,which(substr(names(data_all2),1,1)=="P")]
write.table(data_pbmc,"GSE99254_PBMC.txt",sep="\t",quote=F)

