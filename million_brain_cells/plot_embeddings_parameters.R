source('load_and_save_pcs_tsne.R')

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
gg_color_hue <- function(n) {
	hues = seq(15,375,length=n+1)
	hcl(h=hues, l=65, c=100)[1:n]
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


louvain_gaussian <- as.factor(read.csv('/data/george/Research_Local/million_neurons/louvain_gaussian.csv',header=F)[,1])



set.seed(3)
ridx <- sample.int(length(louvain_gaussian), 5E4)
ridx1E5 <- sample.int(length(louvain_gaussian), 1E5)
ridx2E5 <- sample.int(length(louvain_gaussian), 2E5)

my_theme <- theme( axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
	       axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.title=element_text(size=9, face="bold"), legend.text=element_text(size=9), legend.title=element_text(size=9))

my.cols <- gg_color_hue(length(levels(louvain_gaussian)))
set.seed(1)
my.cols <- sample(my.cols)



#############################################
# Early exaggeration iterations are important
#############################################
toplot <- data.frame( fewExag=iter4k.stop_early_exag_iter250.seed3.vptree[ridx,], 
		     manyExag=iter4k.stop_early_exag_iter2k.seed3.vptree[ridx,],
		     col = as.factor(louvain_gaussian[ridx] ))
g1 <- ggplot(toplot, aes(x = fewExag.1, y = fewExag.2,color=col)) + geom_point(cex=0.1) +guides(color=F) + scale_color_manual(values=my.cols) +my_theme  + ggtitle('Exaggeration for 250 out of 4000 iterations')
g2 <- ggplot(toplot, aes(x = manyExag.1, y = manyExag.2,color=col)) + geom_point(cex=0.1) +guides(color=F)+ scale_color_manual(values=my.cols) + my_theme + ggtitle('Exaggeration for 2000 out of 4000 iterations')
my.cols2 <- my.cols
my.cols2[! levels(louvain_gaussian) %in% c(0,7,14,17,18,19,21,26)] <- "grey"
g3 <- ggplot(toplot, aes(x = fewExag.1, y = fewExag.2,color=col)) + geom_point(cex=0.1) +guides(color=F) + scale_color_manual(values=my.cols2) + my_theme 
g4 <- ggplot(toplot, aes(x = manyExag.1, y = manyExag.2,color=col)) + geom_point(cex=0.1) +guides(color=F)+ scale_color_manual(values=my.cols2) + my_theme 
g <- gridExtra::grid.arrange(g1,g2,g3,g4,nrow=2)
#ggsave(g, file='/data/george/Research_Local/FIt-SNE-paper/figs/early_exag_iter.pdf', width=9,height=9)

ggsave(g, file='/data/george/Research_Local/FIt-SNE-paper/figs/early_exag_iter_raster.jpg', width=180,height=120, unit="mm")

#############################################
# ANNOY works great
#############################################




source('/data/george/Research_Local/FIt-SNE-paper/convenience.R')
toplot <- data.frame( annoy=iter4k.stop_early_exag_iter2k.seed3.annoy[ridx1E5,], 
		     vptree=iter4k.stop_early_exag_iter2k.seed3.vptree[ridx1E5,],
		     col = as.factor(louvain_gaussian[ridx1E5] ))
g1 <- ggplot(toplot, aes(x = annoy.1, y = annoy.2,color=col)) + geom_point(cex=0.1) +guides(color=F) + scale_color_manual(values=my.cols) +my_theme + ggtitle(" ") +get_annotation(toplot, toplot[,1:2],toplot$col,8,3) + ggtitle('Approximate nearest neighbors')
g2 <- ggplot(toplot, aes(x = vptree.1, y = vptree.2,color=col)) + geom_point(cex=0.1) +guides(color=F)+ scale_color_manual(values=my.cols) + my_theme + ggtitle(" ")+get_annotation(toplot, toplot[,3:4],toplot$col,8,3) + ggtitle('Exact nearest neighbors')
g <- gridExtra::grid.arrange(g1,g2,nrow=1)

ggsave(g,file='/data/george/Research_Local/FIt-SNE-paper/figs/vptree_vs_annoy.pdf',width=9, height=4.5)


ggsave(g, file='/data/george/Research_Local/FIt-SNE-paper/figs/vptree_vs_annoy_raster.jpg', width=180,height=90, unit="mm")




#############################################
# We can identify rare populations
#############################################
# Load specific genes
gene_expressions <- list()
fns <- list.files('/data/george/Research_Local/FIt-SNE-paper/million_brain_cells/data/genes')
for (fn in fns) {
	goi <-  strsplit(fn,split='.', fixed=T)[[1]][1]	
	gene_expressions[[goi]] <- read.csv(sprintf("/data/george/Research_Local/FIt-SNE-paper/million_brain_cells/data/genes/%s.csv",goi),header=F)
}


cols <- rep("FALSE",nrow(iter4k.stop_early_exag_iter2k.seed3.annoy))
cols[(gene_expressions[['sncg']] >0) & (gene_expressions[['slc17a8']]>0)] <- "GABAergic subtype (Sncg Slc17a8)    "
cols[(gene_expressions[['spp1']] >0) & (gene_expressions[['col15a1']]>0)] <- "VLMC subtype (Spp1 Col15a1)"
#cols[(gene_expressions[['c1ql2']] >0) & (gene_expressions[['cdh13']]>0) ] <- "L5 PT VISp C1ql2 CDh13"
toplot <- data.frame( iter4k.stop_early_exag_iter2k.seed3.annoy, cols=factor(cols))
g1 <- ggplot(data=toplot[ridx1E5,]) + geom_point( aes(x = X1, y = X2), cex=0.05,alpha=0.5, color="grey") + 
	geom_point(data=subset(toplot,cols!="FALSE"),aes(x = X1, y = X2,color=cols), size=1) +  labs(colour="Cell Type") + labs(x="FIt-SNE 1", y="FIt-SNE 2") + my_theme+ theme(legend.position="none") + ggtitle(' ') 

write.csv(toplot[ridx1E5,1:2],"/data/george/Research_Local/FIt-SNE-paper/figs/Fig1_grey.csv")
write.csv(subset(toplot,cols!="FALSE"),"/data/george/Research_Local/FIt-SNE-paper/figs/Fig1_colors.csv")

toplot <- data.frame( iter4k.stop_early_exag_iter2k.seed3.vptree.50kpts, cols=factor(cols[ridx]))
g3 <- ggplot(toplot, aes(x = X1, y = X2,color=col)) + geom_point( aes(x = X1, y = X2), cex=0.2,alpha=0.5, color="grey") + 
	geom_point(data=subset(toplot,cols!="FALSE"),
	   aes(x = X1, y = X2,color=cols), size=1) +  labs(colour="Cell Type") + 
	guides(colour=guide_legend(override.aes = list(size=5)))+
	labs(x="t-SNE 1", y="t-SNE 2") + my_theme + theme(legend.position="bottom")+ ggtitle(' ')
write.csv(toplot[,1:2],"/data/george/Research_Local/FIt-SNE-paper/figs/Fig1b_grey.csv")
write.csv(subset(toplot,cols!="FALSE"),"/data/george/Research_Local/FIt-SNE-paper/figs/Fig1b_colors.csv")
g<- grid.arrange(g1,g3+theme(legend.position='none'),nrow=1)

ggsave(g,filename="/data/george/Research_Local/FIt-SNE-paper/figs/rare_populations.pdf",width=12,height=6)


legend <- g_legend(g3)
ggsave(legend,filename="/data/george/Research_Local/FIt-SNE-paper/figs/rare_populations_legend.pdf",width=7,height=1)



#############################################
# ANNOY and VP trees of rare population are the same
#############################################

toplot <- data.frame( iter4k.stop_early_exag_iter2k.seed3.annoy, cols=factor(cols))
g1 <- ggplot(data=toplot[ridx1E5,]) + geom_point( aes(x = X1, y = X2), cex=0.05,alpha=0.5, color="grey") + 
	geom_point(data=subset(toplot,cols!="FALSE"),aes(x = X1, y = X2,color=cols), size=1) +  labs(colour="Cell Type") + labs(x="FIt-SNE 1", y="FIt-SNE 2") + my_theme+ theme(legend.position="none")+ggtitle('Approximate nearest neighbors')
toplot <- data.frame( iter4k.stop_early_exag_iter2k.seed3.vptree, cols=factor(cols))

g2 <- ggplot(data=toplot[ridx1E5,]) + geom_point( aes(x = X1, y = X2), cex=0.05,alpha=0.5, color="grey") + geom_point(data=subset(toplot,cols!="FALSE"),aes(x = X1, y = X2,color=cols), size=1) +  labs(colour="Cell Type") + labs(x="FIt-SNE 1", y="FIt-SNE 2") + my_theme+ theme(legend.position="bottom", legend.justification="center")+ ggtitle('Exact nearest neighbors')

g<- grid.arrange(g1,g2+theme(legend.position='none'),nrow=1)
#ggsave(g,filename="/data/george/Research_Local/FIt-SNE-paper/figs/rare_populations_annoy_vp.pdf",width=12,height=6)

legend <- g_legend(g2)
g<- grid.arrange(arrangeGrob( g1,g2+theme(legend.position='none'),nrow=1), legend, nrow=2,heights=c(10,1))

ggsave(g,filename="/data/george/Research_Local/FIt-SNE-paper/figs/rare_populations_annoy_vp_raster.jpg",width=180,height=90, unit="mm")
