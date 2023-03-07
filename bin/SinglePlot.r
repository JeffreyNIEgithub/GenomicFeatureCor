#!/opt/sysoft/R-4.0.2/bin/Rscript
# Usage: Rscript myscript.R
##################################################################################################################
#bed距离画图   包括 相对距离 和 绝对距离（All|overlap|upstream|downstream） 以及 每条染色体上observed和random（overlap）
##################################################################################################################

args<-commandArgs(TRUE)

#################################################
# R包
myPaths <- .libPaths()
new <- c("/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
myPaths <- c(myPaths, new)
.libPaths(myPaths)
library("tidyr", quietly = T,lib.loc="/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
library("ggplot2", quietly = T,lib.loc="/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
library("ggpubr",quietly = T,lib.loc="/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
library("gridExtra",quietly = T,lib.loc="/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
#################################################
# 环境
setwd(args[1])

#################################################
# 数据读入

RelDist=paste(args[2],".reldist.All.tsv",sep="")
relDist <- read.table(RelDist,header=F)
colnames(relDist) <- c("Reldist","Count","Total","Fraction","Type")

AbsDist=paste(args[2],".absdis.tsv",sep="")
absdist <- read.table(AbsDist,header=F)
colnames(absdist) <- c("Chr1","S1","E1","ID1","Chr2","S2","E2","ID2","Compare","Group","Type","Absolute.dist")
absdist$absolute.dist <- (absdist$Absolute.dist)*0.001


AbsDistlog=paste(args[2],".absdis_log.tsv",sep="")
absdistlog <- read.table(AbsDistlog,header=F)
colnames(absdistlog) <- c("Chr1","S1","E1","ID1","Chr2","S2","E2","ID2","Compare","Group","Type","Absolute.dist")
absdistlog$absolute.dist <- log2(absdistlog$Absolute.dist)
#all <- All[,c(1,2,3,4,6)]

# upabsdist <- subset(absdist, absdist$Group == 'upstream')
# downabsdist <- subset(absdist, absdist$Group == 'downstream')
# overlap <- subset(absdist, absdist$Group == 'overlap')
#updownabsdist <- rbind (upabsdist,downabsdist)
#########################
# 简单统计
long <- as.data.frame(table(absdist[which(absdist[c(12)]<0),][c(1,11)]))
widedata <- spread(long,key = "Type",value = "Freq")
overlap=paste(args[2],".overlap.tsv",sep="")
write.table(widedata,overlap,col.names =TRUE,row.names = FALSE,quote=F,sep = "\t") 

cbPalette <- c("#0073C2","#EFC000")
################################################
# 数据可视化
##########################
# 相对距离

# 柱状图
#cbPalette <- c("#0073C2","#EFC000")
pdf1=paste(args[2],".Relative.dist_bar.pdf",sep="")
tltle1=paste(args[2],"-Relative.dist",sep="")
pdf(pdf1,15,10)
ggplot(data=relDist,mapping=aes(x=Reldist,y=Count,fill = Type)) + geom_bar(stat="identity") + facet_grid(Type~.) + labs(title = tltle1,x = NULL) + theme(legend.position="none" ,axis.title = element_text(size=15),axis.text=element_text(size=10),plot.title = element_text(color = 'black', face="bold",size = 32, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 25, angle = 0),axis.text.y = element_text(color = 'black', size = 20, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 25, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text = element_text(color = 'black', size   = 14),strip.text = element_text(face = "bold", size = rel(2.5)))  + scale_fill_manual(values=cbPalette)
dev.off()

##########################
# overlap在各个染色体上的observed及random的count值

# 柱状图
#adPalette <- c("#0073C2","#EFC000")
pdf2=paste(args[2],".overlap_count.dist_bar.pdf",sep="")
tltle2=paste(args[2],"-overlap.count",sep="")
pdf(pdf2,15,10)
ggplot(data=long,mapping=aes(x=Chr1,y=Freq,fill = Type)) + geom_bar(stat="identity") + facet_grid(Type~.) + labs(title = tltle2,x = NULL,y="Count") + theme(legend.position="none" ,axis.title = element_text(size=15),axis.text=element_text(size=10),plot.title = element_text(color = 'black', face="bold",size = 32, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 22, angle = 0),axis.text.y = element_text(color = 'black', size = 20, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 25, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text = element_text(color = 'black', size   = 14),strip.text = element_text(face = "bold", size = rel(2.5)))  + scale_fill_manual(values=cbPalette)
dev.off()

##########################
# 绝对距离

#取log和不去log在一张图
title3=paste(args[2],".Absolute_distance",sep="")
title4=paste(args[2],".Absolute_dist",sep="")
#par(mar = c(15,10, 15,10))
p1 <- ggboxplot(absdist, x="Type", y="absolute.dist", color = "Type", palette = "jco",short.panel.labs = FALSE) + facet_wrap(facets = "Group", scales = "free") + labs(title = title3,x = NULL,y = 'Absolute_distance (Kb)') + theme(plot.title = element_text(color = 'black', face="bold",size = 25, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 18, angle = 0),axis.text.y = element_text(color = 'black', size = 15, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 18, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text  = element_text(color = 'black', size = 14),strip.text = element_text(face = "bold", size = rel(2.0))) + stat_compare_means(method = "t.test",size = 4.2,label.x =1.35) + scale_fill_manual(values=cbPalette)  #按dose进行分面

p2 <- ggboxplot(absdistlog, x="Type", y="absolute.dist", color = "Type", palette = "jco", short.panel.labs = FALSE) + facet_wrap(facets = "Group", scales = "free")  + labs(title = title4,x = NULL,y = 'Log2(Absolute_distance+1) (Kb)') + theme(plot.title = element_text(color = 'black', face="bold",size = 25, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 18, angle = 0),axis.text.y = element_text(color = 'black', size = 18, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 18, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text = element_text(color = 'black', size = 14),strip.text = element_text(face = "bold", size = rel(2.0))) + stat_compare_means(method = "t.test",size = 4.2,label.x =1.35) + scale_fill_manual(values=cbPalette)  #按dose进行分面

pdf3=paste(args[2],".Absolute.dist.pdf",sep="")
pdf(pdf3,width=20,height=10)
grid.arrange(p1,p2,ncol=2)
dev.off()
################################################




