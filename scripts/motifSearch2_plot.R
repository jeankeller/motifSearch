

# R


PATH = "/home/jkeller/work/test_motifSearch/motifSearch_results_20191023-124124/res/"
OUTPATH = "/home/jkeller/work/test_motifSearch/motifSearch_results_20191023-124124/stats/"



require(data.table, quietly=T)
require(tidyr, quietly=T)
library(dplyr, verbose=F)
require(ggplot2, quietly=T)
require(gridExtra, quietly=T)

multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL, heights = NULL, widths = NULL) {
  if (!requireNamespace("gridExtra")) {
    stop("gridExtra package is required. Please install it.")
  }
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots <- length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }
  
  gridExtra::grid.arrange(grobs = plots, newpage = TRUE, layout_matrix = layout, heights = heights, widths = widths)
}

##### Count data #####

df_motif_count <- as.data.frame(fread(paste0(PATH,"motifs_count.txt"), h=T))
motifs <- unique(df_motif_count$motif)


for(i in motifs) {
  plot_list = list()
  subset_motif <- df_motif_count[df_motif_count$motif==i,]
  

  pdf(paste0(OUTPATH,"Hist_motif_plots",i,".pdf"),width = 4, height =8)

  m=matrix(c(1:12), nrow = 6, ncol = 2, byrow = FALSE)
  layout(m,c(1,12),c(9,2,9,2,9,2))
  for(x in c(1:3)){
    par(mar=c(0,0,0,0), ps=8)
    plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE, col="white")
    text(0.5,0.5,"Number of sequences",srt=90,cex=1.5)
    par(mar=c(0,0,0,0), ps=8)
    plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE, col="white")
  }
  
  par(mar=c(0,1.5,1.5,0.5),ps=8,mgp=c(0.6,0.3,0))
  hist(subset_motif$total, breaks=30, main=paste0("Distribution of ",i," in forward and reverse"), xlab="",
                   axes=FALSE,ylab="",cex=1.4, col="gold", xlim=c(0,max(subset_motif$total)))
  abline(v = median(subset_motif$total), col = "black", lwd = 1.2)
  abline(v = quantile(subset_motif$total,0.05), col = "black", lwd = 1,lty=2)
  abline(v = quantile(subset_motif$total,0.95), col = "black", lwd = 1,lty=2)
  legend(x = "topleft",legend =c("Median","Quantiles 5%-95%"), bg="white",col = c("black", "black", "black"),lwd = c(1.2, 1, 1),lty=c(1,2,2), cex=0.8)
  axis(1, las=1)
  axis(2, las=2,tck=-0.01)
  box(lwd=1)
  par(mar=c(0,0,0,0), ps=8)
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE, col="white")
  text(0.5,0.5,"Number of motifs",cex=1.3)
  par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
  hist(subset_motif$count_fwd, breaks=20, main=paste0("Distribution of ",i," in forward"), xlab="",
                   axes=FALSE,ylab="",cex=1.4,col="cornflowerblue", xlim=c(0,max(subset_motif$total)))
  abline(v = median(subset_motif$count_fwd), col = "black", lwd = 1.2)
  abline(v = quantile(subset_motif$count_fwd,0.05), col = "black", lwd = 1,lty=2)
  abline(v = quantile(subset_motif$count_fwd,0.95), col = "black", lwd = 1,lty=2)
  legend(x = "topright",legend =c("Median","Quantiles 5%-95%"), bg="white",col = c("black", "black", "black"),lwd = c(1.2, 1, 1),lty=c(1,2,2), cex=0.8)
  axis(1, las=1)
  axis(2, las=2,tck=0.03)
  box(lwd=1)
  par(mar=c(0,0,0,0), ps=8)
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE, col="white")
  text(0.5,0.5,"Number of motifs",cex=1.3)
  par(mar=c(0,0.7,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
  hist(subset_motif$count_rev, breaks=20, main=paste0("Distribution of ",i," in reverse"), xlab="",
                   axes=FALSE,ylab="",cex=1.4,col="lightcoral", xlim=c(0,max(subset_motif$total)))
  abline(v = median(subset_motif$count_rev), col = "black", lwd = 1.2)
  abline(v = quantile(subset_motif$count_rev,0.05), col = "black", lwd = 1,lty=2)
  abline(v = quantile(subset_motif$count_rev,0.95), col = "black", lwd = 1,lty=2)
  legend(x = "topright",legend =c("Median","Quantiles 5%-95%"), bg="white",col = c("black", "black", "black"),lwd = c(1.2, 1, 1),lty=c(1,2,2), cex=0.8)
  axis(1, las=1)
  axis(2, las=2,tck=0.03)
  box(lwd=1)
  par(mar=c(0,0,0,0), ps=8)
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE, col="white")
  text(0.5,0.5,"Number of motifs",cex=1.3)
  
  dev.off()
  
}


##### Positions data ##### 


fwd <- as.data.frame(fread(paste0(PATH,"motifs_positions_fwd.txt"), h=T))
  fwd$strand <- "fwd"
rev <- as.data.frame(fread(paste0(PATH,"motifs_positions_rev.txt"), h=T))
  rev$strand <- "rev"
fwd_rev <- rbind(fwd,rev)
fwd_rev$mid <- fwd_rev$start+((fwd_rev$stop-fwd_rev$start)/2)

motifs <- unique(fwd_rev$motif)


Size <- 400

window_size <- Size/2
sliding <- seq(from=1, to=max(fwd_rev$stop), by=50)

plot_list = list()





for (motif in motifs) {
  
  subset_motif <- fwd_rev[fwd_rev$motif==motif,]
  plot_list = list()
  YMAX <- vector()
  
  for (j in unique(subset_motif$sequence_id)) {
    
    subset <- subset_motif[subset_motif$sequence_id==j,]
    tot <- vector()

    for (i in sliding){
      win <- subset[subset$mid>i-window_size & subset$mid<i+window_size,]
      density <- nrow(win)
      tot <- append(tot,density)
    }
    
    sliding_tot <- data.frame(x=sliding,y=tot)

    smoothingSpline = smooth.spline(sliding_tot$x,sliding_tot$y, spar=0.07)
    sliding_tot$x <- smoothingSpline$x
    sliding_tot$y <- smoothingSpline$y
    sliding_tot$y <- ifelse(sliding_tot$y<0,0,sliding_tot$y)

    YMAX <- append(YMAX, max(sliding_tot$y))
    rm(sliding_tot)
  }
  
  ylim <- round(max(YMAX)+0.5)
    
  for (j in unique(subset_motif$sequence_id)) {
      
      subset <- subset_motif[subset_motif$sequence_id==j,]
      tot <- vector()
      
      for (slide in sliding){
        win <- subset[subset$mid>slide-window_size & subset$mid<slide+window_size,]
        density <- nrow(win)
        tot <- append(tot,density)
      }
      
      sliding_tot <- data.frame(x=sliding,y=tot)
      
      smoothingSpline = smooth.spline(sliding_tot$x,sliding_tot$y, spar=0.07)
      sliding_tot$x <- smoothingSpline$x
      sliding_tot$y <- smoothingSpline$y
      sliding_tot$y <- ifelse(sliding_tot$y<0,0,sliding_tot$y)
      
      p <- ggplot() + geom_area(sliding_tot, mapping=aes(x=x,y=y)) + scale_y_continuous(limits=c(0,ylim)) +
        labs(title=paste0(j,"_",motif),x="Position in sequence", y = "Total motif number") +
        theme_classic() +
        theme(legend.position="none", plot.margin=margin(t = 10, r = 20, b = 10, l = 20, unit = "pt"),
              plot.title = element_text(size = 10,hjust=0), axis.title = element_text(size = 8))
      
      plot_list[[j]] = p
      
      rm(sliding_tot)
  }
  
  rowplot <- length(unique(subset_motif$sequence_id))
  plotsize <- rowplot/14
  plotsize <- ifelse(plotsize<1,1,round(plotsize+0.5,0))
  
  list <- data.frame(start=1,stop=14)
  
  pdf(paste0(OUTPATH,motif,"_motifs_plots_sliding_windows_",Size,"bp.pdf"),width = 8.27, height = 11.69)
  for (i in 1:plotsize){
    sublist <- list[i,]
    subplot_list <- plot_list[sublist$start:sublist$stop]
    Filter_list <- Filter(Negate(is.null), subplot_list)
    matrixsize <- length(Filter_list)
    layout <- matrix(c(1:14), nrow = 7, ncol = 2, byrow = FALSE)
    multiplot(plotlist = Filter_list, layout = layout)
    list <- rbind(list,data.frame(start=sublist$start+14,stop=sublist$stop+14))
  }
  dev.off()
  
  rm(plot_list)

}



print("Plots -> Done!")

