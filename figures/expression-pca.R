library(componentSkeleton)
rm(list=ls())

ann <- read.delim("/mnt/projects/ikaros/data/samples.csv", stringsAsFactors = F)
ann <- ann[ann$Xeno != "yes" & ann$Exclude != "yes",]
ann$Group[ann$Group=="IK6"] <- "IKN"
ann$Group[ann$Alias=="IKD_5"] <- "IKN"
ann$Group <- factor(ann$Group, levels=c("IKN", "IKD", "IKC"))

#expr <- Matrix.read("/mnt/projects/ikaros/results/anduril/execute/qcReport-samplesClusterHeatmap/vst.csv")
expr <- Matrix.read("/mnt/projects/ikaros/results/anduril/execute/deseqExprMatrix/table.csv")
expr <- expr[,colnames(expr) %in% ann$Alias]
colnames(expr) <- ann$UPN[match(colnames(expr), ann$Alias)]

pdf("/mnt/projects/p2ry8-crlf2/results/figures/expression-pca.pdf")

#for (q in c(0.99, 0.98, 0.97, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0)) {
for (q in c(0)) {

  # filter genes
  #cov <- apply(expr, 1, function(x) (quantile(x, 0.75)-quantile(x, 0.25))/quantile(x, 0.5))
  #expr.filt <- expr[cov >= quantile(cov, q),]
  expr.filt <- expr
  
  # get colors for plotting
  colors <- c('red', 'blue', 'darkgray')
  
  # MDS
  distances <- dist(t(expr.filt))
  MDS <- cmdscale(distances, k=2)
  colnames(MDS) <- c("x","y")
  
  
  # extend axis limits a bit to make sure there is enough space for plotting labels inside
  xspan <- max(MDS[,1])-min(MDS[,1])
  xlim=c(min(MDS[,1])-xspan/10, max(MDS[,1])+xspan/10)
  
  # plot
  par(mar=c(5,2,5,9)) ; par(xpd=TRUE)
  plot(MDS[,1], MDS[,2], type='n', xlab='', ylab='', xlim=xlim, main=sprintf("Top %.0f%% most variable genes (n=%d)", (1-q)*100, nrow(expr.filt)))
  text(MDS[ann$UPN,"x"], MDS[ann$UPN,"y"], ann$UPN, col=colors[as.integer(ann$Group)], cex=0.9)
  legend(x=par("usr")[2]+1, y=max(par("usr")[3:4]), legend=levels(ann$Group), col=colors, pch=19, cex=1.3)
}

dev.off()
