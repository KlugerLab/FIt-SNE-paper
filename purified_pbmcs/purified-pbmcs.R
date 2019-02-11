library(Matrix)
library(cowplot)
library(gridExtra)
library(rsvd)
library(reshape2)
library(ggplot2)

if ( !file.exists('data/purified_pbmc.RData')) {
	urls <- list(b_cells='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz', 
		     cd14_monocytes='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/cd14_monocytes_filtered_gene_bc_matrices.tar.gz',
		     cd34='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd34/cd34_filtered_gene_bc_matrices.tar.gz',
		     cd4_helper='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd4_t_helper/cd4_t_helper_filtered_gene_bc_matrices.tar.gz',
		     regulatory_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/regulatory_t/regulatory_t_filtered_gene_bc_matrices.tar.gz',
		     naive_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_t/naive_t_filtered_gene_bc_matrices.tar.gz',
		     memory_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/memory_t/memory_t_filtered_gene_bc_matrices.tar.gz',
		     cd56_nk='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd56_nk/cd56_nk_filtered_gene_bc_matrices.tar.gz',
		     cytotoxic_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cytotoxic_t/cytotoxic_t_filtered_gene_bc_matrices.tar.gz',
		     naive_cytotoxic='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_cytotoxic/naive_cytotoxic_filtered_gene_bc_matrices.tar.gz'	)


	A <- c()
	labels = c()
	for (i in seq_along(urls)) {
		print(i)
		label <-names(urls)[i]
		fn <- sprintf('data/purified_pbmcs/%s.tar.gz', label)
		download.file(urls[[i]],fn)
		# untar
		fn2 <- sprintf('data/purified_pbmcs/%s_unzipped' ,label)
		untar(fn, exdir=fn2)
		mtx <- as.matrix((readMM(
					 sprintf('%s/filtered_matrices_mex/hg19/matrix.mtx', fn2))))
		genenames <- read.delim(sprintf('%s/filtered_matrices_mex/hg19/genes.tsv', fn2),
					sep = '\t',header = FALSE)[,2]
		rownames(mtx) <- genenames
		if (i>1 && !all(rownames(mtx) == colnames(A))) {
			error('Trying to concatenate a matrix with different genenames')
		}
		A <- rbind(A,t(mtx))
		labels <- c(labels, rep(label,ncol(mtx)))
	}
	save(A,labels,file='data/purified_pbmc.RData')
}else{
	print("Loading");
	load('data/purified_pbmc.RData')
}


num_of_genes_in_cell <- rowSums(A>0)
num_of_cells_in_gene <- colSums(A>0)
keep_cells = which(num_of_genes_in_cell > 400) 
keep_genes = which(num_of_cells_in_gene > 100)

A_ <- A[keep_cells,keep_genes]
labels_ <- as.factor(labels[keep_cells])

exclusive_subset = which(!labels_ %in% c("cd4_helper","cytotoxic_t"))

A_ <- A_[exclusive_subset,]
labels_ <- labels_[exclusive_subset]
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
library(ggplot2)


my_theme <- theme(
axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
plot.title=element_text(size=9, face="bold"), legend.text=element_text(size=9))


g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


labels_ren <- labels_
levels(labels_ren) <- c("B Cells", "CD14 Monocytes", "CD34+", "CD4+ T Helper", "CD56+ NK Cells", "Cytotoxic T Cells", "Memory T Cells", "Naive Cytotoxic T Cells", "Naive T Cells", "Regulatory T Cells")
cbind(levels(labels_), levels(labels_ren))


source('/data/george/Research_Local/FIt-SNE-paper/convenience.R')
df <- data.frame(tsne.seed3.ann,  color=labels_ren)
g1 <- ggplot(df, aes(x = X1, y = X2,color=color)) + geom_point(cex=0.1) + guides(col=guide_legend(nrow=2,override.aes=list(size=2))) + my_theme  +  theme(legend.justification="center",legend.position="bottom",legend.title=element_blank())+ ggtitle("Approximate nearest neighbors") + get_annotation(df, df[,1:2],labels_,8,3)
df <- data.frame(tsne.seed3.vptree,  color=labels_ren)
g2 <- ggplot(df, aes(x = X1, y = X2,color=color)) + geom_point(cex=0.1) +  theme(legend.position="bottom",) + my_theme+ ggtitle("Exact nearest neighbors") + get_annotation(df, df[,1:2], labels_,8,3)
g <- gridExtra::grid.arrange(g1+guides(col=F),g2+guides(col=F),nrow=1)
#ggsave(g,file="/data/george/Research_Local/FIt-SNE-paper/figs/purified_pbmcs.png",width=8,height=4,dpi=600,units="in")

legend <- g_legend(g1)

#ggsave(legend,filename="/data/george/Research_Local/FIt-SNE-paper/figs/purified_pbmcs_legend.pdf",width=9,height=1)

g<- grid.arrange(arrangeGrob( g1+theme(legend.position='none'),g2+theme(legend.position='none'),nrow=1), legend, nrow=2,heights=c(10,2))
ggsave(g,filename="/data/george/Research_Local/FIt-SNE-paper/figs/purified_pbmcs_raster.jpg",width=180,height=90, unit="mm")

