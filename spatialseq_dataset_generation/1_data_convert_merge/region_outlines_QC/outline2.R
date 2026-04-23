#sams <- c("683TL_A2","769BR_D6","769BR_F2","769BR_F4","838BL_A4","B01205A2","B01205A4","B01803A5B5","C01339E1","C01830A5B5","C01934E5F6","C01939E3","C02134A2")
sams <- read.table('sample.15.txt')[,1]
info <- read.table('info2.txt',sep="\t")
donor <- info$V2[match(sams,info$V1)]


library(ggsci)
cols <- unique(c(pal_locuszoom("default")(7),pal_igv("default")(51)))
cols <- cols[-9]

library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)

coll <- read.table('color.type.txt',header=T,sep="\t")

center2center <- 500
diamter <- 220
ratio <- 1e6
binSize <- 1
ruleLen <- 1 #mm
w <- ruleLen*ratio/binSize/center2center/10

addN <- function(x){
        A <- table(x)
        B <- paste(names(A),as.numeric(A),sep=" n=")
        return(B[match(x,names(A))])
}


sums <- array()
labs <- rep(0,nrow(coll))

for(sam in sams){
	su <- readRDS(paste0(sam,".sums.rds"))
	la <- readRDS(paste0(sam,".outline.rds"))
	sums <- rbind(sums,su)
	ll <- table(as.numeric(la))
	ll <- ll[names(ll)!="0"]
	idx <- as.numeric(names(ll))
	labs[idx] <- labs[idx]+ll
}
sums <- sums[-1,]

A <- table(sums$region)
regionCount <- data.frame(region=coll$type,cells=as.numeric(A)[match(coll$type,names(A))],regions=labs)
regionCount$density <- floor(regionCount$cells/regionCount$regions*w*w)


mt <- sums %>% group_by(region) %>% summarise(mean=mean(MT))
mt$mean <- floor(mt$mean+0.5)

prefix <- "total"

pdf(paste0(prefix,".counts.pdf"),width=6,height=6)
#pic <- ggplot(mt,aes(region,mean,fill=region))+geom_bar(stat='identity')+scale_fill_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=mt$region,y=mt$mean,label=paste0(mt$mean,"%"),hjust=0.5,vjust=-1)+ylab("average percentage of MT- genes")+ggtitle(paste(prefix,"- MT gene percentage"))+ylim(c(0,max(mt$mean)*1.05))
#print(pic)
pic <- ggplot(regionCount,aes(region,density,fill=region))+geom_bar(stat='identity')+scale_fill_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=regionCount$region,y=regionCount$density,label=regionCount$density,hjust=0.5,vjust=-1)+ylab(expression("cells per mm"^2))+ggtitle(paste(prefix,"- cell density"))+ylim(c(0,max(regionCount$density)*1.05))
print(pic)
sums$region <- addN(sums$region)
md <- sums %>% group_by(region) %>% summarise(MID_median=median(MIDCount),gene_median=median(geneCount),MT_median=median(MT))

pic <- ggplot(sums,aes(region,MT,color=region))+geom_violin()+geom_boxplot(width=0.2)+scale_colour_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=md$region,y=md$MT_median,label=floor(0.5+md$MT_median),color=cols[1:nrow(md)],hjust=-1.4+0.02*nrow(md),vjust=0.5)+ggtitle(paste(prefix,"- MT gene percentage"))
print(pic)

df <- data.frame(MT=sort(sums$MT))
df$percentage <- seq(nrow(df))*100/nrow(df)
df$range <- floor(df$percentage/5)*5
df$show <- unlist(lapply(seq(nrow(df)),function(x){r <- df$range[x];df$percentage[x]==min(df$percentage[df$range==r])}))
pic <- ggplot(df,aes(percentage,MT))+geom_line()+annotate('point',x=df$percentage[df$show],y=df$MT[df$show])+annotate('text',x=df$percentage[df$show],y=df$MT[df$show]+max(df$MT)*0.025,label=floor(df$MT[df$show]*100)/100,size=3)+theme_bw()+ggtitle(paste(prefix,"- MT gene percentage"))+xlab("cell percentage")+ylab(expression("Mitochondrial gene percentage"))
print(pic)



pic <- ggplot(sums,aes(region,geneCount,color=region))+geom_violin()+geom_boxplot(width=0.2)+scale_colour_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=md$region,y=md$gene_median,label=md$gene_median,color=cols[1:nrow(md)],hjust=-1+0.05*nrow(md),vjust=0.5)+ggtitle(paste(prefix,"- MID counts"))
print(pic)
pic <- ggplot(sums,aes(region,MIDCount,color=region))+geom_violin()+geom_boxplot(width=0.2)+scale_colour_manual(values=cols)+theme_bw()+theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),legend.position = "none")+annotate("text",x=md$region,y=md$MID_median,label=md$MID_median,color=cols[1:nrow(md)],hjust=-1+0.05*nrow(md),vjust=0.5)+ggtitle(paste(prefix,"- gene counts"))
print(pic)
dev.off()
