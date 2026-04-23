##  Rscript ../plotDEG_sigScore_SC.R merged_Border_foamy-Border_ramified.csv
args<-commandArgs(T) 
library(Seurat)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

wd<-"/home/ubuntu/spatial/analysis/MS_brain/_merged_addDon6_3.3/cell_level_regionDEGs/2_combined_donors/"

degfile<- args[1]
deg<-read.csv(paste0(wd, degfile), row.names = "X")
deg<-subset(deg, subset=p_val_adj<0.05)
up_deg<-subset(deg, subset=avg_log2FC>0)
down_deg<-subset(deg, subset=avg_log2FC<0)

up<-rownames(up_deg)
down<-rownames(down_deg)

SeuObj<-readRDS("/home/ubuntu/spatial/analysis/MS_brain/scData/Absinta2021_noCI_noNeurons.rds")
allgenes<-rownames(SeuObj@assays$RNA@data)
up<-intersect(up, allgenes)
down<-intersect(down, allgenes)

SeuObj<-AddModuleScore(SeuObj, features=list(up), name="up")
SeuObj<-AddModuleScore(SeuObj, features=list(down), name="down")

cmpname<-gsub(".csv|merged_", "", degfile)
cmpname<-gsub("Peri-lesion", "Peri_lesion", cmpname)
cmpname<-gsub("_vs_", "-", cmpname)
cmpname<-gsub("_", " ", cmpname)
nameup<-gsub("(.*)[-].*", "\\1", cmpname)
namedown<-gsub(".*[-](.+).*", "\\1", cmpname)

pdf(paste0(gsub(".csv", "", degfile), "_padj005_sigScore_SC.pdf"), width=12, height = 6)
p1<-FeaturePlot(SeuObj, features = paste0("up", "1"), label=T, repel=T, order=T, label.color = "grey", raster=F) +coord_fixed()+
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds"))+
  labs(color="SeuratScore", title = paste0("high in ", nameup))
p2<-FeaturePlot(SeuObj, features = paste0("down", "1"), label=T, repel=T, order=T, label.color = "grey", raster=F) +coord_fixed()+
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "GnBu"))+
  labs(color="SeuratScore", title = paste0("high in ", namedown))
grid.arrange(p1, p2, ncol=2)
dev.off()
