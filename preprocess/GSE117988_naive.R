data_all<-read.delim("GSE117988_raw.expMatrix_PBMC.csv.gz",stringsAsFactors=F)
data_naive<-data_all[,which(!substr(names(data),18,18) %in% c(2,3,4)))]
write.table(data_naive,file="GSE117998_PBMC_naive.txt",quote=F,sep="\t",row.names=F)
