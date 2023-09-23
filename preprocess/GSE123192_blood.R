filelist<-read.delim("GSE123192_filelist.txt",stringsAsFactors=F)
filelist<-filelist[grep("PBMC",filelist[,2]),]

i<-1
	molecules<-read.delim(filelist[i,2],stringsAsFactors=F)
	expression_matrix_df<-as.data.frame(as.matrix(molecules))
	print(dim(molecules))
	print(dim(expression_matrix_df))

for(i in 2:length(filelist[,1])){
	molecules<-read.delim(filelist[i,2],stringsAsFactors=F)
	expression_matrix_df<-merge(expression_matrix_df,as.data.frame(as.matrix(molecules)),by="X")
	print(dim(molecules))
	print(dim(expression_matrix_df))
}
rownames(expression_matrix_df)<-expression_matrix_df[,1]
expression_matrix_df<-expression_matrix_df[,-1]
write.table(expression_matrix_df,file="/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/GEOdata/pretreatment/matrix/GSE123192.txt",quote=F,sep="\t")
