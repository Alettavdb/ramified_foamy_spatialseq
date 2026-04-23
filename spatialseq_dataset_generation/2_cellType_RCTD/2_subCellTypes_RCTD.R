args<-commandArgs(T)

# Rscript /data/work/mydata/runRCTD.R /data/work/mydata/scSub/IMM_sub_clusters.rds /data/work/mydata/spSub/Cellbin_merged_immune-lym.rds
library(Seurat)
library(spacexr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

scData<- args[1] 
spData<- args[2] 

sc.data<-readRDS(scData)
sp.data<-readRDS(spData)
# exclude controls
sp.data<-sp.data[,!grepl("Con", colnames(sp.data))]
sp.data<-JoinLayers(sp.data)

message("RCTD:")
sc.counts<-sc.data@assays$RNA@counts
meta_data<-sc.data@meta.data

# set to default idents
cell_types<-Idents(sc.data); 
names(cell_types)<-rownames(meta_data)
nUMI<-meta_data$nCount_RNA; names(nUMI)<-rownames(meta_data)
reference<-Reference(sc.counts, cell_types, nUMI)

sp.counts<-sp.data@assays$RNA@layers$counts 
rownames(sp.counts)<-rownames(sp.data)
colnames(sp.counts)<-colnames(sp.data)

coords<-sp.data@meta.data[, c("coord_x", "coord_y")]
sp.nUMI<-colSums(sp.counts)
names(sp.nUMI)<-colnames(sp.data)

puck<-SpatialRNA(coords, sp.counts, sp.nUMI)
barcodes <- colnames(puck@counts)
#plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 

myRCTD<-create.RCTD(puck, reference, max_cores = 1, CELL_MIN_INSTANCE = 10, UMI_min = 1,UMI_min_sigma=1) # for parallel processing, the number of cores used. 1 means no parallel processing.
# UMI_min = 1,UMI_min_sigma=1; default only pixels with >100UMI included in analysis
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet") 
saveRDS(myRCTD, file=paste0(spData, "_myRCTD.rds"))

results<-myRCTD@results
table(results$results_df$spot_class)
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA

## combine max weight with meta data
norm_weights = normalize_weights(results$weights); norm_weights<-as.data.frame(norm_weights)
#noReject<-rownames(subset(results$results_df, subset=spot_class!="reject"))
#norm_weights<-norm_weights[which(rownames(norm_weights) %in% noReject),]

maxcol<-colnames(norm_weights)[max.col(norm_weights, ties.method = "first")]
maxweight<-cbind(rownames(norm_weights), maxcol); 
maxweight<-as.data.frame(maxweight)

table(maxweight$maxcol)

sp.data@meta.data$RCTD<-maxweight$maxcol[match(rownames(sp.data@meta.data), maxweight$V1)]
sp.data$RCTD[which(is.na(sp.data$RCTD))]<-"low-quality"
sp.data@meta.data$RCTD_spotclass<-results$results_df$spot_class[match(rownames(sp.data@meta.data), rownames(results$results_df))]
#sp.data <- AddMetaData(sp.data, metadata = myRCTD@results$results_df)
sp.data<-AddMetaData(sp.data, metadata=norm_weights)
sp.data$donor<-substr(rownames(sp.data@meta.data), 1, 4)
saveRDS(sp.data, file=paste0(spData, "_addRCTD.rds"))

write.csv(sp.data@meta.data, file=paste0(spData, "_cell_type_RCTD.meta.csv"))

pdf(paste0(spData, "_cell_type_RCTD.pdf"),width=15,height=6)

regCol<-read.table("/data/work/mydata/region_cols.txt", header=T, sep="\t")
newCol<-regCol[which(regCol$region %in% names(table(sp.data$region))),]
newCol<-newCol[order(newCol$region),]

sp.data@meta.data$region<-factor(sp.data$region, levels=c("Centre ramified", "Border ramified", "Peri-lesion ramified", "NAWM",
                                                "Peri-lesion foamy", "Border foamy", "Centre foamy"))

newCol<-newCol[c(4,2,7,5,6,1,3),]

prctd<-DimPlot(sp.data, reduction="spatial", group.by = "RCTD", split.by="donor",cols = myCol, pt.size = 0.1, ncol=5)+coord_fixed()+theme_void()
grid.arrange(prctd)
p1<-ggplot(sp.data@meta.data,aes(x=RCTD, y=1, fill=RCTD_spotclass))+geom_col(position="fill")+scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = myCol)+labs(x="Sample", y="Proportion", fill="spot_class", title="RCTD spot class")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p2<-ggplot(sp.data@meta.data,aes(x=donor, y=1, fill=RCTD))+geom_col(position="fill")+scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = myCol)+labs(x="Donor", y="Proportion", fill="subtypes", title="Subtypes in donors")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p3<-ggplot(sp.data@meta.data,aes(x=region, y=1, fill=RCTD))+geom_col(position="fill")+scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = myCol)+labs(x="Region", y="Proportion", fill="subtypes", title="Subtypes in regions")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p4<-ggplot(sp.data@meta.data,aes(x=RCTD, y=1, fill=region))+geom_col(position="fill")+scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = newCol$color)+labs(x="Subtypes", y="Proportion", fill="Region", title="Region in subtypes")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p5<-ggplot(sp.data@meta.data,aes(x=RCTD, y=1, fill=region))+geom_col(position="fill")+scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = newCol$color)+labs(x="Subtypes", y="Proportion", fill="Region", title="Region in subtypes")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_grid(.~donor)

p6<-ggplot(sp.data@meta.data,aes(x=region, fill=factor(RCTD)))+geom_bar()+scale_fill_manual(values = myCol)+
    labs(x="Donor", y="Cell count", fill="RCTD maxWeight", title="RCTD subtype counts")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_grid(.~donor)

p7<-ggplot(sp.data@meta.data,aes(x=region, y=1, fill=RCTD))+geom_col(position="fill")+scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = myCol)+labs(x="Region", y="Proportion", fill="subtypes", title="Subtypes in regions")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_grid(.~donor)

grid.arrange(p1, p2, ncol=2)
grid.arrange(p3, p4, ncol=2)
grid.arrange(p5)
grid.arrange(p6)
grid.arrange(p7)
dev.off()
