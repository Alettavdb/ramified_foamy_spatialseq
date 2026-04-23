args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript CellChat_ms_regions.R [DonX] [region] [celltype/all_states]")
}

suppressPackageStartupMessages({
  library(CellChat) # need to be loaded before Seurat
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
})


spData<-"/home/ubuntu/spatial/analysis/MS_brain/_merged_addDon6_3.3/mergedDonors_adj_addStates.rds" # mergedDonors_adj.rds
myDon<-args[1]
myRegion<-args[2]
testGrp<-args[3] # celltype / allsub

dir.create(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp))
myDir<-paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp)

SeuObj<-readRDS(spData)
SeuObj<-subset(SeuObj, subset=donor==myDon & region==myRegion)
gc()

if (testGrp=="allsub") {
 SeuObj<-subset(SeuObj, subset=allsub!="low-quality")
}

SeuObj

#### for spatial
spatial.locs<-SeuObj@meta.data[, c("coord_y", "coord_x")]
spatial.locs$coord_y<-(-1)*spatial.locs$coord_y #flip y coord
#spot.diameter is the theoretical spot size microns; spot is the number of pixels that span the diameter of a theoretical spot size
#scale.factors<-list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres) # ~89.5 for Visium V1

#spot.diameter = cell diameter (based on biological background);   spot: average counts of bin1 in cellbin
# for stereo-seq 10X microscope: 1 pixel=0.5um
# for CellChat 1.6.1
#scale.factors<-list(spot.diameter = 10, spot = 20)

#'ratio' is the conversion factor when converting spatial coordinates from Pixels or other units to Micrometers 
# 'tol' can be the the half value of the minimum center-to-center distance.
scale.factors<-data.frame(ratio=0.5, tol=5)
cellChat<-createCellChat(object=SeuObj, group.by=testGrp, assay = "RNA", datatype = "spatial", coordinates = spatial.locs, spatial.factors = scale.factors) # scale.factors in 1.6.1

cellChat

# L-R database
CellChatDB <- CellChatDB.human
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellChat@DB <- CellChatDB.use

# subsetting is necessary even if using the whole database
cellChat <- subsetData(cellChat)

#future::plan("multiprocess", workers = 4) # do parallel
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

#?computeCommunProb for details
cellChat <- computeCommunProb(cellChat, type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.1,
                              contact.dependent = TRUE, contact.range = 100) #VGP dist ~600um; doesn't change

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellChat <- filterCommunication(cellChat, min.cells = 10)

df.net <- subsetCommunication(cellChat) #slot.name = "netP"
## show all pathways
df.net
write.csv(df.net, file=paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_interactions.csv"))

# calculate probablities
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)


#All the signaling pathways showing significant communications
cellChat@netP$pathways

# Compute the network centrality scores
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") 

# save
saveRDS(cellChat, file = paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_spatial.rds"))

#################
# visualize
groupSize <- as.numeric(table(cellChat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_plots.pdf"), width=15, height=15)
netVisual_circle(cellChat@net$count, vertex.weight = rowSums(cellChat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = rowSums(cellChat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# chord plot of all LR pairs
netVisual_chord_gene(cellChat, sources.use = NULL, targets.use = NULL, lab.cex = 0.8,legend.pos.y = 30)
dev.off()

pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_interaction_heat.pdf"), width=6, height=6)
netVisual_heatmap(cellChat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellChat, measure = "weight", color.heatmap = "Blues")
dev.off()

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellChat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellChat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellChat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellChat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

system("rm *.png")

# pathway heatmaps plots
pheat<-lapply(cellChat@netP$pathways, function(x){
  netVisual_heatmap(cellChat, signaling = x, color.heatmap = "Reds")
  #netVisual_aggregate(cellChat, signaling = x, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "outgoing", vertex.size.max = 5, vertex.label.cex = 3.5)
})

pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_Pathways_heatmap.pdf"))
pheat
dev.off()

## single group plots
pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_singleGroups_signal.pdf"), width=12, height=15)
mat <- cellChat@net$weight
par(mfrow = c(4,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# spatial plots
psigsp<-lapply(cellChat@netP$pathways, function(x){
  netVisual_aggregate(cellChat, signaling = x, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "outgoing", vertex.size.max = 5, vertex.label.cex = 3.5)
})

pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_Pathways_spatial.pdf"))
psigsp
dev.off()

# circle plots
psigcir<-lapply(cellChat@netP$pathways, function(x){
  netVisual_aggregate(cellChat, signaling = x, layout = "circle")
})

pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_CellChat_Pathways_circle.pdf"))
psigcir
dev.off()


pdf(paste0(myDon, "_", gsub(" ", "_", myRegion), "_", testGrp, "_Outgoing-Incoming_patterns.pdf"), width=12, height=6)
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming")
ht1+ht2
dev.off()

pathways.show <- c("CDH")
netVisual_aggregate(cellChat, signaling = pathways.show, layout = "circle")
netAnalysis_signalingRole_network(cellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# bigger circle indicates larger incoming signaling
netVisual_aggregate(cellChat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "outgoing", vertex.size.max = 5, vertex.label.cex = 3.5)

spatialFeaturePlot(cellChat, features = c("APP","CD74"), point.size = 0.8, color.heatmap = "Reds", direction = 1)

