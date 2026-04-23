# neighbors of lymphocytes in MS

allmeta<-read.csv("/home/ubuntu/spatial/analysis/MS_brain/_merged_addDon6_3.3/mergedDonors_adj.rds.meta.csv", row.names = "X")
library(DR.SC)
library(dplyr)

# "Don1A2" "Don2D6" "Don4F2" "Don1F4" "Don7AB" "Don3F4" "Don6AB" "Don8EF" "Don8E3" "Don6A2"

i<-"Don8E3"
  message(i)
  mymeta<-allmeta[allmeta$id==i,]
  typec<-as.data.frame(table(mymeta$celltype))
  totallym<-typec[which(typec$Var1=="lymphocytes"),"Freq"]
  colData<-mymeta[,c("coord_x", "coord_y")]
  colData<-mymeta[,c("coord_x", "coord_y")]
  colnames(colData)<-c("row", "col")
  
  # for stereoseq: 1pixel=0.5um
  nb.test<-getAdj_manual(as.matrix(colData), radius = 100)
  rownames(nb.test)<-rownames(colData)
  colnames(nb.test)<-rownames(colData)
  
  metadata<-mymeta[,c("coord_x", "coord_y", "celltype")]
  metadata$cellnames<-paste0(rownames(metadata), "#", metadata$celltype)
  colnames(nb.test)<-metadata$celltype
  rownames(nb.test)<-metadata$cellnames
  
  message("Lymphocytes neighbors:")
  lym.test<-nb.test[grepl("lymphocytes", rownames(nb.test)),]
  neigh.vas<-lym.test[,grepl("vascular", colnames(lym.test))] # take lymphocytes neighbors
  neigh.vas<-as.data.frame(neigh.vas)
  neigh.lym<-lym.test[,grepl("lymphocytes", colnames(lym.test))] # take vascular neighbors
  neigh.lym<-as.data.frame(neigh.lym)
  
  # sum neighbors
  lymsum<-rowSums(neigh.lym); lymsum<-as.data.frame(lymsum) 
  vassum<-rowSums(neigh.vas); vassum<-as.data.frame(vassum)
  
 neighborsum<-cbind(lymsum, vassum) 
 nonvas.lym<-nrow(neighborsum[which(neighborsum$lymsum<=3 & neighborsum$vassum==0),])
  message(paste0(nonvas.lym, " / ", totallym, "\n"))
  
  neighborsum$lymType<-"vascular"
  neighborsum$lymType[which(neighborsum$lymsum<=3 & neighborsum$vassum==0)]<-"non-vascular"
  rownames(neighborsum)<-gsub("#lymphocytes", "", rownames(neighborsum))
  write.csv(neighborsum, file=paste0(i, "lymType.csv"))
# 19/214
# previous 89/214

# for i in Don1A2 Don1F4 Don2D6 Don3F4 Don4F2 Don6A2 Don6AB Don7AB Don8EF Don8E3; do Rscript check_Neighbors.R $i; done

# cat *.csv --> all donors
