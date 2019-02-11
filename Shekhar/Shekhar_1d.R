source("/data/Linderman/FIt-SNE-paper/Shekhar/class.R")
load("/data/Linderman/FIt-SNE-paper/Shekhar/bipolar_data.Rdata") # https://github.com/broadinstitute/BipolarCell2016
set.seed(3)
print(dim(bipolar_dge))
mt.genes = grep("mt-", rownames(bipolar_dge), value = TRUE)
cells.use = colnames(bipolar_dge)[colSums(bipolar_dge[mt.genes, ])/colSums(bipolar_dge) < 0.1]
bipolar_dge = bipolar_dge[, cells.use]
dim(bipolar_dge)
dsq.bip=scDrop(count.data=bipolar_dge)
dsq.bip=initialize(dsq.bip, min.genes = 500, min.cells = 30, min.counts=60)
dsq.bip@pca.scores = pca.scores
dsq.bip@pca.load = pca.load
data.plot = dsq.bip@pca.scores
data.plot$group = dsq.bip@group
ggplot(data.plot, aes(x = PC1, y = PC2)) + geom_point(aes(colour = factor(group),
                                                          size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
#pca.scores is 27499 by 100, that's the number of cells
#pca.load is 13166 by 100, that's the number of filtered genese
############ 
##  t-SNE Heatmaps
############
source("/data/Linderman/FIt-SNE/fast_tsne.R",chdir=T);
system.time(shekhar_tsne_1d <- fftRtsne(as.matrix(dsq.bip@pca.scores[,1:37]), dims=1,rand_seed = 3))


#We don't need it for the heatmaps, but since we are used to looking at 2D
#t-SNE, let's compute that too 
system.time(shekhar_tsne_2d <- fftRtsne(as.matrix(dsq.bip@pca.scores[,1:37]),  rand_seed = 3))


library(dbscan)
clustering <- dbscan(shekhar_tsne_2d,eps=2, minPts = 40)
ggplot(data.frame(V1=shekhar_tsne_2d[,1], V2=shekhar_tsne_2d[,2], color=as.factor(clustering$cluster)), aes(x=V2,y=V1, color=color) ) +  geom_point(size=0.05)  +
  guides(colour = guide_legend(override.aes = list(size=5), nrow=3)) + theme(legend.position="bottom",title = element_blank())


querygenes <- c("Tacr3","Rcvrn","Syt2","Irx5","Irx6","Vsx1","Hcn4","PkarIIb","Grik1","Gria1","Kcng4", "Hcn1", "Cabp5", "Grm6", "Isl1", "Scgn", "Otx2", "Vsx2","Car8","Sebox","Prkca"); 

fake_genes <- data.frame();
for (i in 0:max(clustering$cluster)){
 fake_gene <- rep(0,ncol(dsq.bip@data))
 fake_gene[clustering$cluster==i] <- 1
 fake_genes <- rbind(fake_genes, t(as.data.frame(fake_gene)) )
}
rownames(fake_genes) <- paste('cluster_', 0:max(clustering$cluster),sep = "")
colnames(fake_genes) <- colnames(dsq.bip@data)
data_with_fake <-rbind(dsq.bip@data, fake_genes)

source('/data/Linderman/t-SNE-Heatmaps_Experimental/t-SNE-Heatmaps/tsnehm_experimental.R')



enrich_with <- 25;
hmtsneout4 <- tsnehm(data_with_fake, c(querygenes, rownames(fake_genes)), shekhar_tsne_1d,clustering$cluster,enrich=enrich_with)
hmtsneout4$heatmap


# Check if all the discovered markers show up in the enriched genes
shekhar_discovered_markers <- c("Pcdh17","Pcdh10", "Erbb4", "Nnat","Col11a1","Sox6","Chrm2","Slitrk5","Lrrtm1","Cck","Lect1","Igfn1","Serpini1","Cpne9","Vstm2b","Casp7")
shekhar_discovered_markers %in% unlist(hmtsneout4$enriched_genes)
all(shekhar_discovered_markers %in% unlist(hmtsneout4$enriched_genes))



#Get the order
hmtsneout4_geneorder <- rev(hmtsneout4$heatmap$x$layout$yaxis2$ticktext)

write.csv(as.matrix(hmtsneout4$matrix[hmtsneout4_geneorder,]),"/data/Linderman/FIt-SNE-paper/Shekhar/Figure2_Data.csv")

goi <- hmtsneout4_geneorder

############################
##  Check with standard hiercarchical clusterp
############################
library(fastcluster)
library(heatmaply)

# Compute distances between cells in PCA space
dd.sub <- dist((as.matrix(dsq.bip@pca.scores[,1:37])))

# Hierarchical clustering on those distances
clust.res <- fastcluster::hclust(dd.sub)

# cutting the tree at 100, which is the number of bins we use in the t-SNE heatmaps (by default)
memb <- cutree(clust.res, k=100)
data.t <- t( rbind( dsq.bip@data[goi, ], fake_genes))
data.agg <- aggregate(data.t, by=list(memb), function(x) mean (x))
my_palette <- colorRampPalette(c("white", "red"))(n = 1000)
#Using the same genes and gene order
toplot <-as.matrix(t(data.agg[,goi]))
toplot <- t(t(toplot) / rowSums(t(toplot)))
slope=100;
intercept=0.01;
toplot2 <- 1/(1+exp(slope*intercept-slope*toplot))
hm.agg.sfig1 <- heatmaply(toplot2,  
           showticklabels = c(FALSE,TRUE),   col=my_palette, dendrogram='column', titleX =FALSE, RowV=FALSE,
          show_legend=FALSE,hide_colorbar = TRUE,fontsize_row = 5,margins = c(70,50,NA,0))
hm.agg.sfig1

## Entirely independent of t-SNE
#library(FNN)
#data_t <- t(dsq.bip@data)
#euclidean_closest <- get.knnx(dsq.bip@data[,], rbind (fake_genes,dsq.bip@data[querygenes[querygenes %in% rownames(dsq.bip@data)],],k=25))
#goi_euclidean <- unique(rownames(dsq.bip@data)[c(euclidean_closest$nn.index)])
#
#data.goi_euclidean.t <- t( dsq.bip@data[goi_euclidean, ])
#data.goi_euclidean.agg <- aggregate(data.goi_euclidean.t, by=list(memb), function(x) mean (x))
#
#toplot <- as.matrix(t(data.goi_euclidean.agg[,-1]))
#toplot <- t(t(toplot) / rowSums(t(toplot)))
#slope=50;
#intercept=0.05;
#toplot2 <- 1/(1+exp(slope*intercept-slope*toplot))
#hm.agg.sfig1_indep <- heatmaply(toplot2,  
#           showticklabels = c(FALSE,TRUE),   col=my_palette, dendrogram='both', titleX =FALSE, RowV=FALSE,
#          show_legend=FALSE,hide_colorbar = TRUE,fontsize_row = 5,margins = c(70,50,NA,0))
#hm.agg.sfig1
#hm.agg.sfig1_indep








## Cutting the tree at 25, which is the number of clusters in the t-SNE
#memb <- cutree(clust.res, k=25)
#data.t <- t( dsq.bip@data[goi, ])
#data.agg <- aggregate(data.t, by=list(memb), function(x) mean (x))
#hm.agg.sfig2 <- heatmaply(as.matrix(t(data.agg[,hmtsneout6_geneorder])),  
#           showticklabels = c(FALSE,TRUE),   col=my_palette, dendrogram='column', titleX =FALSE, RowV=FALSE,
#          show_legend=FALSE,hide_colorbar = TRUE,fontsize_row = 5,margins = c(70,50,NA,0))
#hm.agg.sfig2


