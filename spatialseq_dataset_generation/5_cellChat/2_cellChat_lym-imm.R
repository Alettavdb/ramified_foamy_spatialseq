## lymphocytes cell chat
lym<-readRDS("MS_data/lymphocytes_clustering/lymphocytes_sct20_clustering.rds")
merged<-readRDS("MS_data/mergedDonors_adj_addStates.rds")

# subset immune and lymphocytes
lym_mg<-merged[,merged$celltype=="immune" | merged$celltype=="lymphocytes"]

## add subcluster data
lym$seurat_clusters<-paste0("C", lym$seurat_clusters) # cellchat doesn't allow numbers
lym_mg$lymSub_mg<-as.character(lym$seurat_clusters)[match(rownames(lym_mg@meta.data), rownames(lym@meta.data))]
lym_mg$lymSub_mg[is.na(lym_mg$lymSub_mg)]<-"microglia"

lym_mg$allStates_sub<-as.character(lym$seurat_clusters)[match(rownames(lym_mg@meta.data), rownames(lym@meta.data))]
lym_mg$allStates_sub[is.na(lym_mg$allStates_sub)]<-lym_mg$all_states

table(lym_mg$lymSub_mg)
table(lym_mg$allStates_sub)

saveRDS(lym_mg, file="MS_data/202511_CellChat/immune_lymphocytes-subClusters.rds")

## The rest similar to 1
