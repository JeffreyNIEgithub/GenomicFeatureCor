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
library("ggpubr", quietly = T,lib.loc="/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
library("gridExtra",quietly = T,lib.loc="/home/nieshuai/R/x86_64-pc-linux-gnu-library/4.0")
#################################################
# 环境
setwd(args[1])

#################################################
# 数据读入
RelDist=paste(args[3],"_others.reldist.All.tsv",sep="")
relDist <- read.table(RelDist,header=F)
colnames(relDist) <- c("Reldist","Count","Total","Group","Type")

overlapnum <- read.table("overlap_number",header=T,sep = "\t")

overlapnumlog <- read.table("overlap_number",header=T,sep = "\t")
overlapnumlog$number <- log2(overlapnumlog$number)

overlapboxplot <- read.table("overlap_boxplot_v2.tsv",header=F,sep = "\t")
colnames(overlapboxplot) <- c("Chr1","S1","E1","ID1","Chr2","S2","E2","ID2","Compare","Group","Type","Absolute.dist")
overlapboxplot$absolute.dist <- log2(overlapboxplot$Absolute.dist)

cbPalette <- c("#0073C2","#EFC000")

################################################
# 数据可视化
##########################
# 相对距离

# 柱状图
cbPalette <- c("#0073C2","#EFC000")
pdf=paste(args[3],".Relative.dist_bar.pdf",sep="")
tltle=paste(args[3],"-Relative.dist",sep="")
pdf(pdf,30,20)
ggplot(relDist, aes(x=Reldist, y=Count, fill=Type)) + geom_bar(position=position_dodge(), stat="identity") + facet_wrap(facets = "Group", scales = "free") + labs(title = tltle,x = NULL) + theme(legend.position = "top",axis.title = element_text(size=15),axis.text=element_text(size=10),plot.title = element_text(color = 'black', face="bold",size = 32, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 18, angle = 0),axis.text.y = element_text(color = 'black', size = 20, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 25, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text = element_text(color = 'black', size   = 14),strip.text = element_text(face = "bold", size = rel(1.3))) + scale_fill_manual(values=cbPalette)
dev.off()

#折线图
pdf1=paste(args[3],".Relative.dist_line.pdf",sep="")
tltle1=paste(args[3],"-Relative.dist",sep="")
pdf(pdf1,30,20)
ggplot(relDist, aes(x=Reldist, y=Count, color=Type,shape=Type)) + geom_line(size=0.8) +geom_point(size=2) + facet_wrap(facets = "Group", scales = "free") + labs(title = tltle1,x = NULL) + theme(legend.position = "top",axis.title = element_text(size=15),axis.text=element_text(size=10),plot.title = element_text(color = 'black', face="bold",size = 32, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 18, angle = 0),axis.text.y = element_text(color = 'black', size = 20, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 25, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text = element_text(color = 'black', size   = 14),strip.text = element_text(face = "bold", size = rel(1.3)))  + scale_color_manual(values=cbPalette)
dev.off()
################################################
# 绝对距离

##展示所有repeat（或其他）中谁起到的影响大
#取log和不去log在一张图
title3=paste(args[3],".vs.others Overlap_number",sep="")
title4=paste(args[3],".vs.others Overlap_number",sep="")
#par(mar = c(15,10, 15,10))
p1 <- ggplot(overlapnum, aes(Type, number)) + geom_line(aes(colour = Group, group = Group),size=0.8) + geom_point(size=2, shape=20)+ labs(title = title3,x = NULL) + theme(axis.text.x  = element_text(color = 'black', size = 8, face="bold", angle = 90,hjust = 0.5),axis.text.y = element_text(color = 'black', size = 12, angle = 0),axis.title.y  = element_text(color = 'black', size = 15, face="bold"),legend.position = "top",plot.title = element_text(hjust = 0.5,size=16)) + scale_color_manual(values=cbPalette)
p2 <- ggplot(overlapnumlog, aes(Type, number)) + geom_line(aes(colour = Group, group = Group),size=0.8) + geom_point(size=2, shape=20)+ labs(title = title4,x = NULL,y='Log2(overlap_number)') + theme(axis.text.x  = element_text(color = 'black', size = 8, face="bold", angle = 90,hjust = 0.5),axis.text.y = element_text(color = 'black', size = 12, angle = 0),axis.title.y  = element_text(color = 'black', size = 15, face="bold"),legend.position = "top",plot.title = element_text(hjust = 0.5,size=16)) + scale_color_manual(values=cbPalette)

pdf3=paste(args[3],".overlap_number.pdf",sep="")
pdf(pdf3,width=10,height=15)
grid.arrange(p1,p2,ncol=1)
dev.off()

# pdf("overlap_number_line.pdf",10,8)
# title2=paste(args[3],"-Relative_distance",sep="")
# ggplot(overlapnum, aes(Type, number)) + geom_line(aes(colour = Group, group = Group),size=0.8) + geom_point(size=2, shape=20)+ labs(title = title2,x = NULL) + theme(axis.text.x  = element_text(color = 'black', size = 8, face="bold", angle = 90,hjust = 0.5),axis.text.y = element_text(color = 'black', size = 12, angle = 0),axis.title.y  = element_text(color = 'black', size = 15, face="bold"),legend.position = "top",plot.title = element_text(hjust = 0.5,size=16))
# dev.off()

##overlap箱线图
pdf2=paste(args[3],".overlap_boxplot.pdf",sep="")
title2=paste(args[3],".overlap_boxplot",sep="")
pdf(pdf2,30,20)
ggboxplot(overlapboxplot, x="Type", y="absolute.dist", color = "Type", palette = "jco", short.panel.labs = FALSE) + facet_wrap(facets = "Compare", scales = "free")  + labs(title = title2,x = NULL,y = 'Log2(Absolute_distance+1) (Kb)') + theme(plot.title = element_text(color = 'black', face="bold",size = 25, hjust = 0.5),axis.text.x = element_text(color = 'black', size = 18, angle = 0),axis.text.y = element_text(color = 'black', size = 18, angle = 0),axis.title.x = element_text(color = 'black', size = 18, face="bold",angle = 0),axis.title.y  = element_text(color = 'black', size = 18, face="bold", angle = 90),legend.title = element_text(color = 'black', size  = 14),legend.text = element_text(color = 'black', size = 14),strip.text = element_text(face = "bold", size = rel(1.5))) + stat_compare_means(method = "t.test",size = 4.2,label.x =1.35) + scale_fill_manual(values=cbPalette) #按dose进行分面
dev.off()

################################################


