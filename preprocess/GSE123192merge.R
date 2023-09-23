filelist<-read.delim("filelist.txt",stringsAsFactors=F)
filelist<-filelist[grep("PBMC",filelist[,2]),]

i<-1
	molecules<-read.delim(filelist[i,2],stringsAsFactors=F)
	expression_matrix_df<-as.data.frame(as.matrix(molecules))
	print(dim(expression_matrix_df))

for(i in 2:length(sample)){
	molecules<-read.delim(filelist[i,2],stringsAsFactors=F)
	expression_matrix_df<-cbind(expression_matrix_df,as.data.frame(as.matrix(molecules)))
	print(dim(expression_matrix_df))
}
write.table(expression_matrix_df,"GSE123192_PBMC.txt",sep="\t",quote=F)

