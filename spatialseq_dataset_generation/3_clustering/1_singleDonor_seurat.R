args<-commandArgs(T)

library(Seurat)
#library(ggsci)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(ggpmisc)
library(cowplot)
library(viridis)
library(ggpubr)

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
                brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

wd<-"/data/work/mydata/cellbin/"
myFile<-args[1]
seuData<-paste0(wd, myFile)
SeuObj<-readRDS(seuData)
SeuObj

message("SCT")
SeuObj %>% SCTransform(verbose = FALSE,variable.features.n = 3000, vars.to.regress=c("percent.mt")) %>%
    RunPCA(verbose = FALSE,assay="SCT") %>%
    RunUMAP( dims = 1:20, verbose = FALSE)%>%
    FindNeighbors( dims = 1:20, verbose = FALSE)%>%
    FindClusters(res=0.5, verbose = FALSE) -> SeuObj

saveRDS(SeuObj, file=paste0(seuData, "_sct.rds"))

pdf(paste0(seuData, "_DimPlot.pdf"), height=6, width=9)
p1<-DimPlot(SeuObj, label=T, cols=myCol)+coord_fixed()+NoLegend()
p2<-DimPlot(SeuObj, reduction="spatial", cols=myCol)+theme_void()+coord_fixed()
#p2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=seurat_clusters))+geom_tile()+theme_void()+coord_fixed()+scale_fill_manual(values = myCol)

grid.arrange(p1, p2, ncol=2)
dev.off()  

# markers
cluster_markers_all <- FindAllMarkers(object = SeuObj, assay = "RNA", only.pos = TRUE)
write.csv(cluster_markers_all,paste0(seuData,"_deg.csv"))

myRedu<-"spatial"
myGrp<-"seurat_clusters"

pdf(paste0(seuData, "_singleClusters_", myRedu, ".pdf"), height=4, width=5)
m_data <- SeuObj@meta.data
clst <- as.data.frame(table(m_data[, myGrp]))
m_data<-cbind(m_data, SeuObj$umap@cell.embeddings)
min_x <- min(m_data$coord_x)
max_x <- max(m_data$coord_x)
min_y <- min(m_data$coord_y)
max_y <- max(m_data$coord_y)

for (i in clst$Var1){
  nClst<-clst$Freq[clst$Var1==i]
  n_data <- m_data[which(m_data$seurat_clusters != i),]
  n_data$seurat_clusters = as.character(n_data$seurat_clusters)
  n_data$seurat_clusters = "background"
  n_data$seurat_clusters = as.factor(n_data$seurat_clusters)
  o_data <- m_data[which(m_data$seurat_clusters == i),]
  if (myRedu=="umap") {
    p1 <- ggplot()+
    geom_point(data = n_data,aes(x = as.numeric(umap_1), y = as.numeric(umap_2)),size = 0.3, color = "grey90")+
    geom_point(data = o_data,aes(x = as.numeric(umap_1), y = as.numeric(umap_2)),size = 0.3, color = myCol[as.numeric(i)+1])+
    theme_classic()+
    coord_fixed()+
    ggtitle(paste0("cluster ",i, " (", nClst, ")"))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    scale_x_continuous(limits = c(min_x, max_x))+
    scale_y_continuous(limits = c(min_y, max_y))
  } else {
    p1 <- ggplot()+
    geom_point(data = n_data,aes(x = as.numeric(coord_x), y = as.numeric(coord_y)),size = 0.3, color = "grey90")+
    geom_point(data = o_data,aes(x = as.numeric(coord_x), y = as.numeric(coord_y)),size = 0.3, color = myCol[as.numeric(i)+1])+
    theme_classic()+
    coord_fixed()+
    ggtitle(paste0("cluster ",i, " (", nClst, ")"))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    scale_x_continuous(limits = c(min_x, max_x))+
    scale_y_continuous(limits = c(min_y, max_y))
  }

  df <- data.frame("Gene"=head(cluster_markers_all$gene[cluster_markers_all$cluster==i],n=20),
                   "padj"=format(head(cluster_markers_all$p_val_adj[cluster_markers_all$cluster==i],n=20),digits=2))
  if (as.numeric(count(df)) == 0){
    Gene <- "NA"
    padj <- "NA"
    df = data.frame(Gene,pValue)
    p2 <- ggplot()+
      annotate("table",size=3,x=0,y=0,label=list(df))+
      theme_void()+
      theme(plot.title = element_text(size = 20, face = "bold"))
  }else{
    p2 <- ggplot()+
      annotate("table",size=3,x=0,y=0,label=list(df))+
      theme_void()+
      theme(plot.title = element_text(size = 20, face = "bold"))
  }
  p.list <- list(p1,p2)
  print(plot_grid(plotlist = p.list,ncol = 2,rel_widths = c(4,2)))
}
dev.off()

