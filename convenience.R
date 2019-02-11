require(FNN)
 get_annotation <- function(df, tsne,truelabels,xshift = 5, textsize=4) {
	nn.ann <- get.knn(tsne,k=1)
	nn.ann.lab <- truelabels[nn.ann$nn.index]
	m1 <- mean(nn.ann.lab == truelabels)
	print(100*m1)
	annotate(x=max(df[,1])-xshift,
		y=min(df[,1]), 
		label=sprintf("1N Error = %.1f%%", (1-m1)*100),
		size=textsize,geom="text")
 }
