
library(Matrix)
GSE148190<-read.delim("filelist.txt",stringsAsFactors=F)
GSE148190<-GSE148190[grep("blood",GSE148190[,2]),]

sample<-unique(substr(GSE148190[,2],1,10))
print(length(sample))
i<-1
        group<-GSE148190[grep(sample[i],GSE148190[,2]),2]
        cellBarcode<-read.table(paste0("GSE148190/",group[grep("barcodes",group)]))
        geneName<-read.table(paste0("GSE148190/",group[grep("genes",group)]))
        molecules<-Matrix::readMM(paste0("GSE148190/",group[grep("matrix",group)]))
        rownames(molecules) <- geneName[,1]
        colnames(molecules) <- paste0("MM_",i,"_",cellBarcode[,1])
        geneMat<-as.matrix(molecules)
        print(dim(geneMat))

for(i in 2:length(sample)){
        group<-GSE148190[grep(sample[i],GSE148190[,2]),2]
        cellBarcode<-read.table(paste0("GSE148190/",group[grep("barcodes",group)]))
        geneName<-read.table(paste0("GSE148190/",group[grep("genes",group)]))
        molecules<-Matrix::readMM(paste0("GSE148190/",group[grep("matrix",group)]))
        rownames(molecules)<-geneName[,1]
        colnames(molecules)<-paste0("MM_",i,"_",cellBarcode[,1])
        geneMat<-cbind(geneMat,as.matrix(molecules))
        print(dim(geneMat))
}

geneInfo<-read.delim("Human_GRCh38.87.Gene_info.txt")
row.names(geneMat)<-geneInfo[match(row.names(geneMat),geneInfo[,1]),2]

library(limma)
geneMat2<-geneMat[-which(is.na(rownames(geneMat))),]
geneMat2<-avereps(geneMat2)
geneMat2[1:4,1:4]
dim(geneMat2)
write.table(geneMat2,"GSE148190_PBMC.txt",sep="\t",quote=F)


