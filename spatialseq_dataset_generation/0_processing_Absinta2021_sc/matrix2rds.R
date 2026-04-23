
# Rscript /data/work/MS_brain/scData/matrix2rds.R 
setwd("/data/work/MS_brain/scData")
library(data.table)
library(Seurat)
counts<-fread("GSE180759_expression_matrix.csv.gz")
rownames(counts)<-counts$V1 ## automatically added colname by fread
counts$V1<-NULL

mtx<-as.matrix(counts)
rownames(mtx)<-rownames(counts)

meta<-read.delim("GSE180759_annotation.txt", sep="\t", row.names="nucleus_barcode")
#meta<-subset(meta, subset=rownames(meta) %in% colnames(counts))

SeuObj<-CreateSeuratObject(counts=mtx, meta.data=meta)
SeuObj

saveRDS(SeuObj, file="GSE180759_seu.rds")

# remove chronic inactive and neurons
SeuObj<-SeuObj[,!grepl("chronic_inactive", SeuObj$pathology)]
SeuObj<-subset(SeuObj, subset=cell_type!="neurons")
message("removing inactive and neurons:")
SeuObj

SeuObj<-NormalizeData(SeuObj)
saveRDS(SeuObj, file="Absinta2021_noCI_noNeurons.rds")
