
library(Matrix)
GSE155698<-read.delim("filelist.txt",stringsAsFactors=F)
GSE155698<-GSE155698[grep("PBMC",GSE155698[,2]),]
sample<-unique(substr(GSE155698[,2],1,10))
print(length(sample))

sample<-dir("GSE155698")
sample<-sample[grep("PDAC",sample)]
i<-1
	cellbarcodes <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/barcodes.tsv.gz"))
	genenames <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/features.tsv.gz"))
	molecules <- Matrix::readMM(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/matrix.mtx.gz"))
	rownames(molecules) <- genenames[,1]
	colnames(molecules) <- paste0("PDAC_",i,"_",cellbarcodes[,1])
	expression_matrix_df <-as.matrix(molecules)
	print(dim(expression_matrix_df))

for(i in 2:length(sample)){
        cellbarcodes <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/barcodes.tsv.gz"))
        genenames <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/features.tsv.gz"))
        molecules <- Matrix::readMM(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/matrix.mtx.gz"))
	rownames(molecules) <- genenames[,1]
	colnames(molecules) <- paste0("PDAC_",i,"_",cellbarcodes[,1])
	expression_matrix_df <- cbind(expression_matrix_df,as.matrix(molecules))
	print(dim(expression_matrix_df))
}

geneInfo<-read.delim("Human_GRCh38.87.Gene_info.txt.gz")
geneMat<-as.matrix(expression_matrix_df)
row.names(geneMat)<-geneInfo[match(row.names(geneMat),geneInfo[,1]),2]

library(limma)
geneMat2<-geneMat[-which(is.na(rownames(geneMat))),]
geneMat2<-avereps(geneMat2)
geneMat2[1:4,1:4]
dim(geneMat2)
write.table(geneMat2,"GSE155698_PDAC.txt",sep="\t",quote=F)








library(Matrix)
setwd("/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/GEOdata/pretreatment/matrix")
GSE155698<-read.delim("/zfssz3/pub/database/ftp.ncbi.nih.gov/geo/series/GSE155nnn/GSE155698/suppl/filelist.txt",stringsAsFactors=F)
GSE155698<-GSE155698[grep("PBMC",GSE155698[,2]),]
sample<-unique(substr(GSE155698[,2],1,10))
print(length(sample))

sample<-dir("/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/GEOdata/GEOdata/GSE155698")
sample<-sample[grep("Healthy",sample)]
i<-1
        cellbarcodes <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/barcodes.tsv.gz"))
        genenames <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/features.tsv.gz"))
        molecules <- Matrix::readMM(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/matrix.mtx.gz"))
        rownames(molecules) <- genenames[,1]
        colnames(molecules) <- paste0("normal3_",i,"_",cellbarcodes[,1])
        #expression_matrix_df 是该样品的表达矩阵,行是gene，列是cell_barcode
        expression_matrix_df <-as.matrix(molecules)
        print(dim(expression_matrix_df))

for(i in 2:length(sample)){
        cellbarcodes <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/barcodes.tsv.gz"))
        genenames <- read.table(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/features.tsv.gz"))
        molecules <- Matrix::readMM(paste0("../../GEOdata/GSE155698/",sample[i],"/filtered_feature_bc_matrix/matrix.mtx.gz"))
        rownames(molecules) <- genenames[,1]
        colnames(molecules) <- paste0("normal3_",i,"_",cellbarcodes[,1])
        expression_matrix_df <- cbind(expression_matrix_df,as.matrix(molecules))
        print(dim(expression_matrix_df))
}

geneInfo<-read.delim("/hwfssz5/ST_PRECISION/OCG/zhouxiaorui/LC_CTC/pipeline2019/Human_GRCh38.87.Gene_info.txt.gz")
geneMat<-as.matrix(expression_matrix_df)
row.names(geneMat)<-geneInfo[match(row.names(geneMat),geneInfo[,1]),2]

library(limma)
geneMat2<-geneMat[-which(is.na(rownames(geneMat))),]
geneMat2<-avereps(geneMat2)
geneMat2[1:4,1:4]
dim(geneMat2)
write.table(geneMat2,"GSE155698_normal3.txt",sep="\t",quote=F)

