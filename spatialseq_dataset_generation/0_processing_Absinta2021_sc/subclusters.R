SeuObj<-readRDS("MS_brain/scData/GSE180759_seu.rds")

# IMM = immune + lymphocytes
# OL = oligodendrocytes + opc
# OPC = opc

myname<-"opc"

sub<-subset(SeuObj, subset=cell_type==myname) # | cell_type=="opc"
sub
submeta<-read.table("MS_brain/scData/meta_subclustering_OPC_Fig1i.tsv", header=T, sep="\t")
dim(submeta)
table(submeta$seurat_clusters2)

sub$subClusters<-submeta$seurat_clusters2[match(submeta$barcodes, rownames(sub@meta.data))]
table(sub$subClusters)

Idents(sub)<-"subClusters"
saveRDS(sub, file="MS_brain/scData/OPC_sub_clusters.rds")
