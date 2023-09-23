library(Matrix)
GSE128066<-read.delim("filelist.txt",stringsAsFactors=F)
GSE128066<-GSE128066[grep("PBL",GSE128066[,2]),]

sampleID<-unique(substr(GSE128066[,2],1,10))
print(length(sampleID))
i<-1
        group<-GSE128066[grep(sampleID[i],GSE128066[,2]),2]
        cellBarcode<-read.table(group[grep("barcodes",group)])
        geneName<-read.table(group[grep("genes",group)])
        molecules<-Matrix::readMM(group[grep("matrix",group)])
        rownames(molecules) <- geneName[,1]
        colnames(molecules) <- paste0("normal2_",i,"_",cellBarcode[,1])
        #geneMat<-as.data.frame(as.matrix(molecules))
        geneMat<-as.matrix(molecules)
        print(dim(geneMat))

for(i in 2:length(sampleID)){
        group<-GSE128066[grep(sampleID[i],GSE128066[,2]),2]
        cellBarcode<-read.table(group[grep("barcodes",group)])
        geneName<-read.table(group[grep("genes",group)])
        molecules<-Matrix::readMM(group[grep("matrix",group)])
        rownames(molecules)<-geneName[,1]
        colnames(molecules)<-paste0("normal2_",i,"_",cellBarcode[,1])
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
write.table(geneMat2,"GSE128066_PBMC.txt",sep="\t",quote=F)


