# The analysis of Hrvatin et al. (2018) data was inspired by Huang et al.
# (2018), the code for which is available at
# https://github.com/mohuangx/SAVER-paper, and data available at
# https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0
library(gridExtra)
library(ggplot2)
library(Seurat)
library(ggplot2)
library(cowplot)

setwd('/data/george/Research_Local/FIt-SNE-paper/Hrvatin')
if (!file.exists('data/hrvatin_full.rds')) {
	x <- as.matrix(read.csv("data/GSE102827_merged_all_raw.csv", header = TRUE, row.names = 1, check.names = FALSE))
	print(dim(x))
	# Filter out genes with expression less than 0.00003
	x1 <- x[rowMeans(x) >= 0.00003, ]

	# Filter out genes with non-zero expression in less than 4 cells
	x2 <- x1[rowSums(x1 != 0) >= 4, ]
	print(dim(x2))
	x <- x2
	dim(x2)
	saveRDS(x, "data/hrvatin_full.rds")
}else{
	print("Reading")
	x <- readRDS('data/hrvatin_full.rds')
}




print(dim(x))
ct.dat <- read.csv("data/hrvatin_celltypes.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
ct.dat2 <- ct.dat[!is.na(ct.dat$celltype), ]
ct.dat2$celltype[!is.na(ct.dat2$subtype)] <- ct.dat2$subtype[!is.na(ct.dat2$subtype)]
table(ct.dat2$celltype)
#ct.dat3 <- ct.dat2[ct.dat2$celltype != "ExcL23", ]
ct.dat3 <- ct.dat2
ident <- ct.dat3[colnames(x), 4]
x.sub <- x[, !is.na(ident)]
print(dim(x.sub))
ident2 <- ct.dat3[colnames(x.sub), 4]





A_ <- t(x.sub)
totalUMIPerCell <- rowSums(A_);
A_norm <- sweep(A_, 1, totalUMIPerCell, '/');
A_norm <- A_norm * 10E3
A_norm <- log(A_norm +1);

#center the columns
A_rm <- colMeans(A_norm)
A_norm_c <- sweep(A_norm, 2,A_rm) 

library(rsvd)
set.seed(3)
fastDecomp <- rsvd(A_norm_c, 25,q=4);
PCs<- fastDecomp$u %*% diag(fastDecomp$d);

source('/data/george/Research_Local/FIt-SNE_Experimental/linqiaozhi/FIt-SNE/fast_tsne.R', chdir=T)

tsne.seed3.ann <- fftRtsne(as.matrix(PCs), rand_seed = 3, ann_not_vptree=TRUE)
tsne.seed3.vptree <- fftRtsne(as.matrix(PCs), rand_seed = 3,ann_not_vptree=FALSE)


my_theme <- theme(
axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
legend.title=element_blank(), plot.title=element_text(size=9, face="bold"), legend.text=element_text(size=7))

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


source('/data/george/Research_Local/FIt-SNE-paper/convenience.R')



df <- data.frame(tsne.seed3.ann,  color=ident2)
g1 <- ggplot(df, aes(x = X1, y = X2,color=color)) + geom_point(cex=0.1) + guides(col=guide_legend(nrow=4,override.aes=list(size=2))) + my_theme  +    theme(legend.position="bottom",legend.justification="center") + ggtitle("Approximate nearest neighbors") +get_annotation(df, df[,1:2],ident2,10,3)
df <- data.frame(tsne.seed3.vptree,  color=ident2)
g2 <- ggplot(df, aes(x = X1, y = X2,color=color)) + geom_point(cex=0.1) +  theme(legend.position="bottom",) + my_theme + ggtitle("Exact nearest neighbors") + get_annotation(df, df[,1:2],ident2,10,3)
g <- gridExtra::grid.arrange(g1+guides(col=F),g2+guides(col=F),nrow=1)
#ggsave(g,file="/data/george/Research_Local/FIt-SNE-paper/figs/hrvatin.png",width=8,height=4,dpi=600,units="in")



legend <- g_legend(g1)
#ggsave(legend,filename="/data/george/Research_Local/FIt-SNE-paper/figs/hrvatin_legend.pdf",width=12,height=2)


g<- grid.arrange(arrangeGrob( g1+theme(legend.position='none'),g2+theme(legend.position='none'),nrow=1), legend, nrow=2,heights=c(10,3))
ggsave(g,filename="/data/george/Research_Local/FIt-SNE-paper/figs/hrvatin_raster.jpg",width=180,height=90, unit="mm")


