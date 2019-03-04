setwd("~/Documents/GitHub/RNAPoser/manuscript/")

# pose classification 
plot_clscores <- function(file, pred = TRUE, ylab = "", plotname = ""){
  # read in data
  data <- read.table(file, col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))
  data$status <- 0
  
  # chose whether to plot classification scores for the actual best or predicted best
  if(pred){
    data$status[data$cs_pred>0.5] <- 1
    data$cs <- data$cs_pred
    data$rmsd <- data$rmsd_pred
  } else {
      data$status[data$cs_true>0.5] <- 1
      data$cs <- data$cs_true
      data$rmsd <- data$rmsd_true
  }
  
  # setup of color scheme
  data$color <- "red"
  data$color[data$status == 1] <- "forestgreen"
  
  # make plot (note: suppressed automatic tick labels)
  plot(data$rmsd, data$cs, lwd = "1", bg = data$color, pch = 21, xaxt = "n", yaxt = "n", ylim = c(0.0  , 1.0), xlim =c (0, 10), ylab = ylab, xlab = "RMSD (Å)", main = plotname)
  # customize tick labels for x-axis
  # see: https://www.statmethods.net/advgraphs/axes.html
  axis(1, at=seq(0, 10, 2),labels=c("0.0", "2.0", "4.0", "6.0", "8.0", "10.0"), col.axis="black", las=2)
  # customize tick labels for y-axis
  axis(2, at=c(0, 0.25, 0.50, 0.75, 1.0), labels=c("0.00", "0.25", "0.50", "0.75", "1.00"), col.axis="black", las=2)
  abline(h=0.5, lwd = "1", lty = "dotted")
  abline(v=2.5, lwd = "1", lty = "dotted")
  #abline(v=mean(data$rmsd[data$status==0]), lwd = "1", lty = "dashed", col = "red")
  #abline(v=mean(data$rmsd[data$status==1]), lwd = "1", lty = "dashed", col = "forestgreen")
  #points(data$rmsd, data$cs, lwd = "0", bg = data$color, pch = 21, ylim = c(0.0, 1.0), xlim =c (0, 10), ylab = ylab, xlab = "RMSD (Å)", main = plotname)
  return(data)
}

# pose classification 
plot_distribution_clscores <- function(pred = TRUE, ylab = ""){
  d1 <- read.table("poses_classification_scores_xray_loo.txt", col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))
  d2 <- read.table("poses_classification_scores_xray_valid.txt", col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))

  d1$grp <- "nmr"
  d2$grp <- "xray"
  d3 <- rbind(d1, d2)
  d3$grp <- "all"
  data <- rbind(d1, d2, d3)
  
  if(pred){
    data$cs <- data$cs_pred
  } else {
    data$cs <- data$cs_true
  }
  boxplot(cs~grp, data = data, yaxt = "n", ylim = c(0.0, 1.0), col=c("white", "gray7 0", "gray30"), ylab = ylab)
  axis(2, at=c(0, 0.25, 0.50, 0.75, 1.0), labels=c("0.00", "0.25", "0.50", "0.75", "1.00"), col.axis="black", las=2)
  #abline(h=0.5, lwd = "1", lty = "dotted")
  return(data)
}



# new figures

# dev size: 9.777778 10.375000 
dev.new(width=9.777778, height=10.375000)
scale <- 1
shift <- 0.0
shifty <- 0.3
par(fig=(scale)*c(0.0+shift ,0.4+shift, 0.0+shifty, 0.5+shifty), new=FALSE)
nmr <- plot_clscores(file = "poses_classification_scores_xray_loo.txt", pred = TRUE, ylab = "CLscore", plotname = "")
d1 <- read.table("poses_classification_scores_xray_loo.txt", col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))
dy <- 0.30
ya <- 0.265
yb <- ya+dy
par(fig=(scale)*c(0.0+shift ,0.4+shift, ya+shifty, yb+shifty), new=TRUE)
boxplot(d1$rmsd_pred, axes=FALSE, xaxt = "n", ylim = c(0.0, 10.0), ylab = "", las = 2, horizontal = TRUE, outline = FALSE)

par(fig=(scale)*c(0.23+shift, 0.47+shift, 0.0+shifty, 0.5+shifty), new=TRUE)
boxplot(d1$cs_pred, axes=FALSE, yaxt = "n", ylim = c(0.0, 1.0), ylab = "", outline = FALSE)
par(fig=(scale)*c(0.45+shift, 0.85+shift, 0.0+shifty, 0.5+shifty), new=TRUE)
xray <- plot_clscores(file = "poses_classification_scores_xray_valid.txt", pred = TRUE, ylab = "CLscore", plotname = "")
d1 <- read.table("poses_classification_scores_xray_valid.txt", col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))
dy <- 0.30
ya <- 0.265
yb <- ya+dy
par(fig=(scale)*c(0.45+shift, 0.85+shift, ya+shifty, yb+shifty), new=TRUE)
boxplot(d1$rmsd_pred, axes=FALSE, xaxt = "n", ylim = c(0.0, 10.0), ylab = "", las = 2, horizontal = TRUE)
par(fig=(scale)*c(0.23+0.45+shift, 0.85+0.07+shift, 0.0+shifty, 0.5+shifty), new=TRUE)
boxplot(d1$cs_pred, axes=FALSE, yaxt = "n", ylim = c(0.0, 1.0), ylab = "")
dev.copy2pdf(file = "clscore.pdf")
stop()


# plot atomTypes frequency
atm <-read.table("atomTypes.txt", col.names = c("name", "counts"))
atm <- atm[order(atm$counts, decreasing = TRUE), ]
barplot(atm$counts, names.arg = c(as.character(atm$name)), las = 2, horiz = TRUE, ylim = c(25,0), cex.axis = 0.8, cex.names = 0.8, ylab = "Atom Types", xlab = "Counts")
