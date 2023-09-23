
library(Matrix)
GSE139324<-read.delim("/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/GEOdata/rawData/GSE139324_filelist.txt",stringsAsFactors=F)
GSE139324<-GSE139324[grep("PBMC",GSE139324[,2]),]
GSE139324<-GSE139324[grep("HNSCC",GSE139324[,2]),]
sample<-unique(substr(GSE139324[,2],1,10))
print(length(sample))
i<-1
	group<-GSE139324[grep(sample[i],GSE139324[,2]),2]
	cellbarcodes <- read.table(group[grep("barcodes",group)])
	genenames <- read.table(group[grep("genes",group)])
	molecules <- Matrix::readMM(group[grep("matrix",group)])
	rownames(molecules) <- genenames[,1]
	colnames(molecules) <- paste0("HNSCC_",i,"_",cellbarcodes[,1])
	expression_matrix_df <-as.matrix(molecules)
	print(dim(expression_matrix_df))

for(i in 2:length(sample)){
	group<-GSE139324[grep(sample[i],GSE139324[,2]),2]
	cellbarcodes <- read.table(group[grep("barcodes",group)])
	genenames <- read.table(group[grep("genes",group)])
	molecules <- Matrix::readMM(group[grep("matrix",group)])
	rownames(molecules) <- genenames[,1]
	colnames(molecules) <- paste0("HNSCC_",i,"_",cellbarcodes[,1])
	expression_matrix_df <- cbind(expression_matrix_df,as.matrix(molecules))
	print(dim(expression_matrix_df))
}


geneInfo<-read.delim("/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/LC_CTC/pipeline2019/Human_GRCh38.87.Gene_info.txt")
geneMat<-as.matrix(expression_matrix_df)
row.names(geneMat)<-geneInfo[match(row.names(geneMat),geneInfo[,1]),2]

library(limma)
geneMat2<-geneMat[-which(is.na(rownames(geneMat))),]
geneMat2<-avereps(geneMat2)
geneMat2[1:4,1:4]
dim(geneMat2)
write.table(geneMat2,file="GSE139324_PBMC.txt",sep="\t",quote=F)



library(Matrix)
GSE139324<-read.delim("GSE139324_filelist.txt",stringsAsFactors=F)
GSE139324<-GSE139324[grep("PBMC",GSE139324[,2]),]
GSE139324<-GSE139324[grep("HD",GSE139324[,2]),]
sample<-unique(substr(GSE139324[,2],1,10))
print(length(sample))
i<-1
        group<-GSE139324[grep(sample[i],GSE139324[,2]),2]
        cellbarcodes <- read.table(group[grep("barcodes",group)])
        genenames <- read.table(group[grep("genes",group)])
        molecules <- Matrix::readMM(group[grep("matrix",group)])
        rownames(molecules) <- genenames[,1]
        colnames(molecules) <- paste0("normal1_",i,"_",cellbarcodes[,1])
        expression_matrix_df <-as.matrix(molecules)
        print(dim(expression_matrix_df))

for(i in 2:length(sample)){
        group<-GSE139324[grep(sample[i],GSE139324[,2]),2]
        cellbarcodes <- read.table(group[grep("barcodes",group)])
        genenames <- read.table(group[grep("genes",group)])
        molecules <- Matrix::readMM(group[grep("matrix",group)])
        rownames(molecules) <- genenames[,1]
        colnames(molecules) <- paste0("normal1_",i,"_",cellbarcodes[,1])
        expression_matrix_df <- cbind(expression_matrix_df,as.matrix(molecules))
        print(dim(expression_matrix_df))
}


geneInfo<-read.delim("Human_GRCh38.87.Gene_info.txt")
geneMat<-as.matrix(expression_matrix_df)
row.names(geneMat)<-geneInfo[match(row.names(geneMat),geneInfo[,1]),2]

library(limma)
geneMat2<-geneMat[-which(is.na(rownames(geneMat))),]
geneMat2<-avereps(geneMat2)
geneMat2[1:4,1:4]
dim(geneMat2)
write.table(geneMat2,"GSE139324_normal1.txt",sep="\t",quote=F)


