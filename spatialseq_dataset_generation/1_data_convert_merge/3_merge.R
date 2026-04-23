# run with Seurat V5 -- result in merged unfiltered data
# Rscript /data/work/MS_brain/merge15.R Cellbin/Bin50
args<-commandArgs(T) 

library(Seurat)
library(ggplot2)

myType<-args[1]

wd<-"/data/work/MS_brain/"
# read all rds files
if (myType=="Cellbin") {
  dataFiles <- lapply(Sys.glob("/data/work/MS_brain/cellbin/*cellTypes.rds"), readRDS)
} else {
  dataFiles <- lapply(Sys.glob("/data/work/MS_brain/bin50/*regions.rds"), readRDS)
}
dataFiles

sampleInfo<-read.table("/data/work/MS_brain/sample-info.txt", header=T, sep="\t")
merged<-merge(dataFiles[[1]], y=dataFiles[2:15], add.cell.ids=sampleInfo$id) ###!!!

merged$id<-substr(rownames(merged@meta.data),1,6)

# add coord info as embedding
sp.embedding<-merged@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
merged$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')

#merged<-NormalizeData(merged)
saveRDS(merged, file=paste0(wd, myType,"_merged15.rds"))
write.csv(merged@meta.data, file=paste0(wd, myType, "_merged15.meta.csv"))

regCol<-read.table("/data/work/MS_brain/region_cols.txt", header=T, sep="\t")
newCol<-regCol[which(regCol$region %in% names(table(merged$region))),]
newCol<-newCol[order(newCol$region),]

pdf(paste0(wd, myType, "_merged_regions.pdf"))
DimPlot(SeuObj, reduction="spatial", group.by="region", split.by="id", ncol=4, cols=newCol$color, raster=T)+coord_fixed()+theme_void()
dev.off()

#FeaturePlot(merged, features="SPP1", reduction="spatial", split.by="id", raster=F, order=T)+coord_fixed()
