args<-commandArgs(T)

if (length(args)<2) {
  cat("Usage: Rscript singleDEG.R [region 1] [region 2]\n")
  q()
}

library(Seurat)
library(dplyr)
library(ggplot2)

x<-"/home/ubuntu/spatial/analysis/MS_brain/_merged_addDon6_3.3/mergedDonors.rds"

id1<-args[1] #"Border foamy"
id2<-args[2] #"Border ramified"


SeuObj<-readRDS(x)
SeuObj<-SeuObj[,!is.na(SeuObj$region)]
SeuObj<-subset(SeuObj, subset=region!="Grey matter")
SeuObj<-NormalizeData(SeuObj)  
Idents(SeuObj)<-"region"
allReg<-list(names(table(SeuObj$region)))
deg<-FindMarkers(SeuObj, ident.1 = id1, ident.2 = id2)
deg$threshold<-ifelse(deg$avg_log2FC>=0.378 & deg$p_val_adj<0.01,"UP", ifelse(deg$avg_log2FC<=-0.378 & deg$p_val_adj<0.01, "DOWN", "NO"))
write.csv(deg, file=paste0("merged_", gsub(" ", "_", id1), "-", gsub(" ", "_", id2), ".csv"))

message(paste0(id1, " vs ", id2))
table(deg$threshold)

deg$delabel <- rownames(deg)
deg$delabel[abs(deg$avg_log2FC)<1.5 | deg$p_val_adj>0.01]<-NA

pdf(paste0("merged_", gsub(" ", "_", id1), "-", gsub(" ", "_", id2), "_vol.pdf"))
ggplot(as.data.frame(deg), aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold, label=delabel))+geom_point()+
  scale_colour_manual(values=c("UP"="red","DOWN"="blue", "NO"="grey"))+labs(title=paste0(id1, " vs ", id2))+
  geom_text(size=3, nudge_y = 1.8)
dev.off()
