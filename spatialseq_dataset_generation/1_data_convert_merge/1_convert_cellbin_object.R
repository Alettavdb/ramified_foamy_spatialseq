# construct seuobj from Lianghan's separation files
args<-commandArgs(T)

if (length(args)<1) {
  cat("Usage: Rscript convert.R [name]\n")
  q()
}

suppressPackageStartupMessages({
  library(Seurat)
  library(RColorBrewer)
  library(gridExtra)
  library(ggplot2)
  library(ggpubr)
})

myCol<-unique(c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"),brewer.pal(12, "Paired"), 
                brewer.pal(8, "Accent"), brewer.pal(8, "Dark2")))

lhwd<-"seuobj/"
rctdwd<-"intermediate/"
name<-args[1]

message(name)

spData<-paste0(lhwd, name, ".cb.seuobj.rds")
spData2<-paste0("MS_big/", gsub("_", "", name), ".cb.seuobj.rds")
infotxt<-paste0(rctdwd, name, ".RCTD.diff.celltype.region.txt")

if (file.exists(spData)){
  SeuObj<-readRDS(spData)
} else {
  SeuObj<-readRDS(spData2)
}


SeuObj
DefaultAssay(SeuObj)<-"RNA"
SeuObj[['SCT']]<-NULL

SeuObj<-NormalizeData(SeuObj)

#message("original metadata:")
#head(SeuObj@meta.data)

mymeta<-read.table(infotxt, header=T, sep="\t")  
## correspond x, y to SeuObj cell coordinates ####
mymeta$x<-mymeta$x+min(SeuObj$coor_x)
mymeta$y<-mymeta$y+min(SeuObj$coor_y)
rownames(mymeta)<-paste0(mymeta$x, "_", mymeta$y)
mymeta$x<-NULL
mymeta$y<-NULL

message("Cell type and region counts: ")
table(mymeta$region)
table(mymeta$celltype)

# combine metadata
SeuObj<-AddMetaData(SeuObj, mymeta)
colnames(SeuObj@meta.data)[which(colnames(SeuObj@meta.data)=="coor_x")]<-"coord_x"
colnames(SeuObj@meta.data)[which(colnames(SeuObj@meta.data)=="coor_y")]<-"coord_y"

SeuObj$orig.ident<-name

#message("added metadata:")
#head(SeuObj@meta.data)

# add coord info as embedding
sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')

saveRDS(SeuObj, file=paste0(name, "_cellTypes.rds"))

message("Converted cell type and region counts: ")
table(SeuObj$region)
table(SeuObj$celltype)

#p0.1<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=nFeature_RNA))+geom_tile()+theme_void()+coord_fixed()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "RdPu")) 
p0.1<-FeaturePlot(SeuObj, reduction="spatial", features="nFeature_RNA", raster=T)+theme_void()+coord_fixed()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"))
#p0.2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=nCount_RNA))+geom_tile()+theme_void()+coord_fixed()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "RdPu")) 
p0.2<-FeaturePlot(SeuObj, reduction="spatial", features="nCount_RNA", raster=T)+theme_void()+coord_fixed()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"))

pvio1<-ggviolin(SeuObj@meta.data, x="orig.ident", y="nFeature_RNA", add = "boxplot", fill="skyblue")+ggtitle("nGene")
pvio2<-ggviolin(SeuObj@meta.data, x="orig.ident", y="nCount_RNA", add = "boxplot", fill="skyblue4")+ggtitle("nCount")
pvio3<-ggviolin(SeuObj@meta.data, x="orig.ident", y="percent.mt", add = "boxplot", fill="yellow")+ggtitle("percent.mt")

# match region color to defined
regCol<-read.table("region_cols.txt", header=T, sep="\t")
newCol<-regCol[which(regCol$region %in% names(table(SeuObj$region))),]
newCol<-newCol[order(newCol$region),]

p1<-DimPlot(SeuObj, reduction="spatial", group.by="celltype", cols=myCol, raster=F)+theme_void()+coord_fixed()
p2<-DimPlot(SeuObj, reduction="spatial", group.by="region", cols=newCol$color, raster=F)+theme_void()+coord_fixed()

p1.2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=celltype))+geom_tile()+theme_void()+coord_fixed()+scale_fill_manual(values = myCol)
p2.2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=region))+geom_tile()+theme_void()+coord_fixed()+scale_fill_manual(values = newCol$color) 

pdf(paste0(name, "_cellTypes_dim.pdf"), width=9, height=4)
grid.arrange(p0.1, p0.2, ncol=2, top="Spatial distribution of nGene and nCount")
grid.arrange(pvio1, pvio2, pvio3, ncol=3)
grid.arrange(p1, p2, ncol=2, top=name)
grid.arrange(p1.2, p2.2, ncol=2, top=name)
dev.off()
