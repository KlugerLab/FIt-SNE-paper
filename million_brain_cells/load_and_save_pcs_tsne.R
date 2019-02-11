library(rsvd) 
library(ggplot2)
source('/data/george/Research_Local/FIt-SNE_Experimental/linqiaozhi/FIt-SNE/fast_tsne.R', chdir=T)

# Run scanpy_preprocessing.py first
if ( !file.exists( 'data/processed_neurons.RData')) {
	A_norm <- read.csv('data/processed_neurons.csv',header=F)
	save(A_norm,file='data/processed_neurons.RData')
}else{
	load('data/processed_neurons.RData')
}
#load('data/processed_neurons.RData')


#############
# PCA and t-SNE
#############
if ( !file.exists( 'data/processed_neurons_PCs.RDS')) {
	    A_rm <- colMeans(A_norm)
	    A_norm_c <- sweep(A_norm, 2,A_rm) 
	    system.time ( fastDecomp <- rsvd(A_norm_c, 100,q=4) )  #2157 seconds
	    colSums(A_norm_c[,1:5])
	    rm(A_norm_c)
	    gc()
	    PCs<- fastDecomp$u %*% diag(fastDecomp$d);
	    saveRDS(PCs, file='data/processed_neurons_PCs.RDS' )
}else{
	PCs <- readRDS(file='data/processed_neurons_PCs.RDS')
}



tsne_fn <- '/data/george/Research_Local/FIt-SNE-paper/million_brain_cells/data/iter4k.stop_early_exag_iter2k.seed3.vptree.RDS'
if ( !file.exists( tsne_fn)) {
	time.results <- system.time ( iter4k.stop_early_exag_iter2k.seed3.vptree <- fftRtsne(as.matrix(PCs[,1:50]),
											max_iter=4000,rand_seed = 3,
											stop_early_exag_iter = 2000,
											ann_not_vptree=FALSE)) 
	saveRDS(iter4k.stop_early_exag_iter2k.seed3.vptree, file=tsne_fn)
} else{
	iter4k.stop_early_exag_iter2k.seed3.vptree <- readRDS(tsne_fn)
}

tsne_fn <- '/data/george/Research_Local/FIt-SNE-paper/million_brain_cells/data/iter4k.stop_early_exag_iter250.seed3.vptree.RDS'
if ( !file.exists( tsne_fn)) {
	time.results <- system.time ( iter4k.stop_early_exag_iter250.seed3.vptree <- fftRtsne(as.matrix(PCs[,1:50]),
											 max_iter=4000,rand_seed = 3,
											 stop_early_exag_iter = 250,
											 ann_not_vptree=FALSE)) 
	saveRDS(iter4k.stop_early_exag_iter250.seed3.vptree, file=tsne_fn)
} else{
	iter4k.stop_early_exag_iter250.seed3.vptree <- readRDS(tsne_fn)
}


tsne_fn <- '/data/george/Research_Local/FIt-SNE-paper/million_brain_cells/data/iter4k.stop_early_exag_iter2k.seed3.annoy.RDS'
if ( !file.exists( tsne_fn)) {
	system.time ( iter4k.stop_early_exag_iter2k.seed3.annoy <- fftRtsne(as.matrix(PCs[,1:50]),
								max_iter=4000,rand_seed = 3,
								stop_early_exag_iter=2000,
								ann_not_vptree=TRUE)
	)
	saveRDS(iter4k.stop_early_exag_iter2k.seed3.annoy, file = tsne_fn)
} else{
	iter4k.stop_early_exag_iter2k.seed3.annoy <- readRDS(tsne_fn)
}


#fofo
tsne_fn <- '/data/george/Research_Local/FIt-SNE-paper/million_brain_cells/data/iter4k.stop_early_exag_iter2k.seed3.vptree.50kpts.RDS'
if ( !file.exists( tsne_fn)) {
	set.seed(3)
	ridx <- sample.int(nrow(PCs), 5E4)
	system.time ( iter4k.stop_early_exag_iter2k.seed3.vptree.50kpts <- fftRtsne(as.matrix(PCs[ridx,1:50]),
								max_iter=4000,rand_seed = 3,
								stop_early_exag_iter=2000,
								ann_not_vptree=FALSE)
	)
	saveRDS(iter4k.stop_early_exag_iter2k.seed3.vptree.50kpts, file = tsne_fn)
} else{
	iter4k.stop_early_exag_iter2k.seed3.vptree.50kpts <- readRDS(tsne_fn)
}






#system.time ( bh_test <- fftRtsne(as.matrix(PCs[,1:50]), max_iter=max_iter,rand_seed = 3, stop_early_exag_iter = 2000, ann_not_vptree=FALSE, fft_not_bh=F,nthreads=1)) 
