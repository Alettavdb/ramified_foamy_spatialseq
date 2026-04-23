## subclustering lym
Idents(SeuObj)<-"celltype"

lym<-subset(SeuObj, subset=celltype=="lymphocytes")
lym<-NormalizeData(lym)
lymtype<-read.csv("MS_data/allCombined_lymType.csv", row.names="X")
lym$lymType<-lymtype$lymType[match(rownames(lym@meta.data), rownames(lymtype))]

saveRDS(lym, file="MS_data/lymphocytes_rna.rds")

## sct dim20
lym %>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>%
  RunPCA(verbose = FALSE,assay="SCT") %>%
  RunUMAP( dims = 1:20, verbose = FALSE)%>%
  FindNeighbors(dims = 1:20, verbose = FALSE)%>%
  FindClusters(verbose = F, res=0.5) -> lym

lym<-FindClusters(lym, res=0.2)

DefaultAssay(lym) <- "RNA"
saveRDS(lym, file="MS_data/lymphocytes_sct20_clustering.rds")

seuData<-"MS_data/lymphocytes_sct"

pdf(paste0(seuData, "20_DimPlot.pdf"), height=6, width=9)
p1<-DimPlot(lym, label=T, cols=myCol)+coord_fixed()+labs(title="Seurat clusters")
p2<-DimPlot(lym, reduction="spatial", cols=myCol)+theme_void()+coord_fixed()
p3<-DimPlot(lym, group.by = "region")+coord_fixed()
p4<-DimPlot(lym, group.by="lymType")+coord_fixed()
#p2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=seurat_clusters))+geom_tile()+theme_void()+coord_fixed()+s                                                       cale_fill_manual(values = myCol)
grid.arrange(p1)
grid.arrange(p2)
grid.arrange(p3)
grid.arrange(p4)
dev.off()

# markers
cluster_markers_all <- FindAllMarkers(object = lym, assay = "RNA", only.pos = TRUE)
write.csv(cluster_markers_all,paste0(seuData,"20res0.2_deg.csv"))


## re clustering dim 10
DefaultAssay(lym)<-"SCT"
lym %>%
  RunUMAP( dims = 1:10, verbose = FALSE)%>%
  FindNeighbors(dims = 1:10, verbose = FALSE)%>%
  FindClusters(verbose = F, res=0.5) -> lym

DimPlot(lym)+coord_fixed()
DimPlot(lym, group.by = "all_states")+coord_fixed()

DefaultAssay(lym)<-"RNA"
lym<-NormalizeData(lym)
FeaturePlot(lym, features = "CD4", coord.fixed = T)

pdf(paste0(seuData, "10_DimPlot.pdf"), height=6, width=9)
p1<-DimPlot(lym, label=T, cols=myCol)+coord_fixed()+labs(title="Seurat clusters")
p2<-DimPlot(lym, reduction="spatial", cols=myCol)+theme_void()+coord_fixed()
p3<-DimPlot(lym, group.by = "region")+coord_fixed()
p4<-DimPlot(lym, group.by="lymType")+coord_fixed()
#p2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x),as.numeric(coord_y),fill=seurat_clusters))+geom_tile()+theme_void()+coord_fixed()+s                                                       cale_fill_manual(values = myCol)
grid.arrange(p1)
grid.arrange(p2)
grid.arrange(p3)
grid.arrange(p4)
dev.off()

# markers
cluster_markers_all <- FindAllMarkers(object = lym, assay = "RNA", only.pos = TRUE)
write.csv(cluster_markers_all,paste0(seuData,"10res0.5_deg.csv"))



###########################

# single clusters
myRedu<-"umap"
myGrp<-"seurat_clusters"

pdf(paste0(seuData, "_singleClusters_", myRedu, ".pdf"), height=4, width=5)
m_data <- lym@meta.data
clst <- as.data.frame(table(m_data[, myGrp]))
m_data<-cbind(m_data, lym$umap@cell.embeddings)
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
