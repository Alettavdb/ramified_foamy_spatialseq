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

lhwd<-"intermediate/"
seuwd<-"seuobj_bin50/"
name<-args[1]

message(name)

spData<-paste0(seuwd, name, ".bin50.seuobj.rds")
SeuObj<-readRDS(spData)
SeuObj
DefaultAssay(SeuObj)<-"RNA"
#SeuObj[['SCT']]<-NULL

SeuObj<-NormalizeData(SeuObj)

lab <- readRDS(paste0(lhwd,name,".outline.rds"))
ctype <- read.table('script/color.type.txt',header=T,sep="\t")
ctype <- c(NA,ctype$type)

s <- dim(lab)
xs <- floor(SeuObj$coor_x*s[2]/max(SeuObj$coor_x))
xs <- xs-min(xs)+1
ys <- floor(SeuObj$coor_y*s[1]/max(SeuObj$coor_y))
ys <- ys-min(ys)+1

xs[xs>dim(lab)[2]] <- dim(lab)[2]
ys[ys>dim(lab)[1]] <- dim(lab)[1]
SeuObj$region <- ctype[1+unlist(lapply(seq(xs),function(x){lab[ys[x],xs[x]]}))]

colnames(SeuObj@meta.data)[which(colnames(SeuObj@meta.data)=="coor_x")]<-"coord_x"
colnames(SeuObj@meta.data)[which(colnames(SeuObj@meta.data)=="coor_y")]<-"coord_y"

SeuObj$orig.ident<-name

# add coord info as embedding
sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')

saveRDS(SeuObj, file=paste0(name, "_BIN50_regions.rds"))

head(SeuObj@meta.data)

message("Converted cell type and region counts: ")
table(SeuObj$region)

p0.1<-FeaturePlot(SeuObj, reduction="spatial", features="nFeature_RNA", raster=T)+theme_void()+coord_fixed()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"))
p0.2<-FeaturePlot(SeuObj, reduction="spatial", features="nCount_RNA", raster=T)+theme_void()+coord_fixed()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "Reds"))

pvio1<-ggviolin(SeuObj@meta.data, x="orig.ident", y="nFeature_RNA", add = "boxplot", fill="skyblue")+ggtitle("nGene")
pvio2<-ggviolin(SeuObj@meta.data, x="orig.ident", y="nCount_RNA", add = "boxplot", fill="skyblue4")+ggtitle("nCount")
pvio3<-ggviolin(SeuObj@meta.data, x="orig.ident", y="percent.mt", add = "boxplot", fill="yellow")+ggtitle("percent.mt")

# match region color to defined
regCol<-read.table("region_cols.txt", header=T, sep="\t")
newCol<-regCol[which(regCol$region %in% names(table(SeuObj$region))),]
newCol<-newCol[order(newCol$region),]

p2.2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=region))+geom_tile()+theme_void()+coord_fixed()+scale_fill_manual(values = newCol$color) 

pdf(paste0(name, "_BIN50_dim.pdf"), width=9, height=4)
grid.arrange(p0.1, p0.2, ncol=2, top="Spatial distribution of nGene and nCount")
grid.arrange(pvio1, pvio2, pvio3, ncol=3)
grid.arrange(p2.2, ncol=2, top=name)
dev.off()
