args <- commandArgs();
if(length(args)!=9){
        write(paste("Rscript",unlist(strsplit(args[4],"="))[2],"<FI:png>","<FI:seuobj>","<INT:binSize>","<STR:prefix>",sep=" "),stderr());
        q();
}

filePNG <- args[6]
fileSeu <- args[7]
binSize <- args[8]
prefix <- args[9]
#filePNG <- '769BR_D6.RNA.1_10.png';fileSeu <-"bin50/769BR_D6.bin50.seuobj.rds";binSize <- 50;prefix <- "769BR_D6.bin50"
#filePNG <- '683TL_A2.RNA.1_10.m.png';fileSeu <-"683TL_A2.cb.seuobj.rds";binSize <- 1;prefix <- "683TL_A2"

ratio_png <- 10
ratio_diff <- 1

center2center <- 500
diamter <- 220
ratio <- 1e6
binSize <- 1
ruleLen <- 1 #mm
w <- ruleLen*ratio/binSize/center2center/10


info <- read.table('info2.txt',sep="\t")
donor <- prefix
if(donor %in% info$V1){donor <- info$V2[info$V1==donor]}

#library(pheatmap)
library(png)
library(ggplot2)
loc <- readPNG(filePNG)
col <- read.table('color.type.txt',header=T,sep='\t')
col$R <- col$R/255
col$G <- col$G/255
col$B <- col$B/255



s <- dim(loc)
if(T){
dif <- matrix(0.5,nrow=s[1],ncol=s[2])
lab <- matrix(0,nrow=s[1],ncol=s[2])

for(n in seq(nrow(col))){
	write(col[n,4],file=stderr())
for(i in seq(s[1])){
	for(j in seq(s[2])){
		scr <- 3
		if(loc[i,j,4]>0){
			scr <- (abs(loc[i,j,1]-col[n,1])+abs(loc[i,j,2]-col[n,2])+abs(loc[i,j,3]-col[n,3]))
		}
		if(scr<dif[i,j]){
			dif[i,j] <- scr
			lab[i,j] <- n
		}
	}
}
	
}
#pheatmap(dif,cluster_rows=F,cluster_cols=F)

saveRDS(lab,paste(prefix,"outline.rds",sep="."))
}
lab <- readRDS(paste(prefix,"outline.rds",sep="."))

library(ggsci)
cols <- unique(c(pal_locuszoom("default")(7),pal_igv("default")(51)))
cols <- cols[-9]

library(Seurat)
library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)


seuobj <- readRDS(fileSeu)
seuobj <-  subset(seuobj, features =row.names(seuobj@assays$RNA)[-1*grep("^mt|^Mt-|^MT-",row.names(seuobj@assays$RNA))])
xs <- floor(seuobj$coor_x*s[2]/max(seuobj$coor_x))
xs <- xs-min(xs)+1
ys <- floor(seuobj$coor_y*s[1]/max(seuobj$coor_y))
ys <- ys-min(ys)+1

xs[xs>dim(lab)[2]] <- dim(lab)[2]
ys[ys>dim(lab)[1]] <- dim(lab)[1]
seuobj$region <- unlist(lapply(seq(xs),function(x){lab[ys[x],xs[x]]}))

addN <- function(x){
	A <- table(x)
	B <- paste(names(A),as.numeric(A),sep=" n=")
	return(B[match(x,names(A))])
}

sums <- array()
for(rg in seq(nrow(col))){
	idx <- which(seuobj$region==rg)
	if(length(idx)>0){
		expr <- seuobj@assays$RNA[,idx]
		sums <- rbind(sums,data.frame(region=rg,MIDCount=apply(expr,2,sum),geneCount=apply(expr,2,function(x){sum(x>0)}),MT=seuobj$percent.mt[idx]))
	}
}
sums <- sums[-1,]
sums$region <- col$type[sums$region]
sums$sample <- prefix
saveRDS(sums,paste0(prefix,".sums.rds"))
lab_density <- table(as.numeric(lab))
names(lab_density) <- unlist(lapply(as.numeric(names(lab_density)),function(x){if(x==0){"Others"}else{col$type[x]}}))
regionCount <- table(sums$region)
regionCount <- floor(regionCount/lab_density[match(names(regionCount),names(lab_density))]*w*w)
regionCount <- data.frame(region=names(regionCount),density=as.numeric(regionCount))

mt <- data.frame(region=seuobj$region,percentage=seuobj$percent.mt)
mt <- mt %>% group_by(region) %>% summarise(mean=mean(percentage))
mt$mean <- floor(mt$mean+0.5)
mt <- mt[mt$region>0,]
mt$region <- col$type[mt$region]


pdf(paste0(prefix,".counts.pdf"),width=6,height=6)
#pic <- ggplot(mt,aes(region,mean,fill=region))+geom_bar(stat='identity')+scale_fill_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=mt$region,y=mt$mean,label=paste0(mt$mean,"%"),hjust=0.5,vjust=-1)+ylab("average percentage of MT- genes")+ggtitle(paste(donor,"- MT gene percentage"))+ylim(c(0,max(mt$mean)*1.05))
#print(pic)
pic <- ggplot(regionCount,aes(region,density,fill=region))+geom_bar(stat='identity')+scale_fill_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=regionCount$region,y=regionCount$density,label=regionCount$density,hjust=0.5,vjust=-1)+ylab(expression("cells per mm"^2))+ggtitle(paste(donor,"- cell density"))+ylim(c(0,max(regionCount$density)*1.05))
print(pic)
sums$region <- addN(sums$region)
#md <- sums %>% group_by(region) %>% summarise(MID_median=median(MIDCount),gene_median=median(geneCount))

md <- sums %>% group_by(region) %>% summarise(MID_median=median(MIDCount),gene_median=median(geneCount),MT_median=median(MT))

pic <- ggplot(sums,aes(region,MT,color=region))+geom_violin()+geom_boxplot(width=0.2)+scale_colour_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=md$region,y=md$MT_median,label=floor(0.5+md$MT_median),color=cols[1:nrow(md)],hjust=-1.4+0.02*nrow(md),vjust=0.5)+ggtitle(paste(prefix,"- MT gene percentage"))
print(pic)


pic <- ggplot(sums,aes(region,geneCount,color=region))+geom_violin()+geom_boxplot(width=0.2)+scale_colour_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=md$region,y=md$gene_median,label=md$gene_median,color=cols[1:nrow(md)],hjust=-1+0.05*nrow(md),vjust=0.5)+ggtitle(paste(donor,"- MID counts"))
print(pic)
pic <- ggplot(sums,aes(region,MIDCount,color=region))+geom_violin()+geom_boxplot(width=0.2)+scale_colour_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=md$region,y=md$MID_median,label=md$MID_median,color=cols[1:nrow(md)],hjust=-1+0.05*nrow(md),vjust=0.5)+ggtitle(paste(donor,"- gene counts"))
print(pic)
dev.off()
