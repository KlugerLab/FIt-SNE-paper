totalits = 1000;
bh_errors <- rep(0, totalits)
fft_errors <- rep(0, totalits)
for (it in 0:(totalits-1)){
  bh_gradient <- read.csv(sprintf('temp/bh_gradient%d.txt',it), header=F)[,-1]
  fft_gradient <- read.csv(sprintf('temp/fft_gradient%d.txt',it),header=F)[,-1]
  exact_gradient <- read.csv(sprintf('temp/exact_gradient%d.txt',it), header=F)[,-1]
   bh_errors[it+1] <- mean(colMeans((abs(bh_gradient - exact_gradient)/abs(exact_gradient))))
   fft_errors[it+1] <- mean(colMeans((abs(fft_gradient - exact_gradient)/abs(exact_gradient))))
}

require(cowplot)
theme_set(theme_minimal())
noerror <- which(0==bh_errors)
bh_errors[noerror] <- NA;
toplot <- data.frame(it = 1:totalits, errors <- (c(bh_errors, fft_errors)), group = c(rep("BH", totalits), rep("FI", totalits)))
acc <- ggplot(toplot, aes(x = it, y = errors)) + geom_point(aes(colour = factor(group)), size=0.5) +
  labs(x="Iteration", y="Mean Relative Error")  +scale_y_log10(breaks=c( 1e-12, 1e-9,  1e-7, 1e-5, 1e-3, 1e-1)) + theme(legend.title = element_blank()) 
#9 rows removed of the first 30
acc

#ggsave("../../fast_tsne_applications_paper/nature_methods_revision2/accuracy2.pdf",width = 7.5, height=5,units = "in",acc) #5x 3
ggsave("../../fast_tsne_applications_paper/nature_methods_revision2/accuracy2.jpg",width = 18, height=9,units = "cm",acc) #5x 3

