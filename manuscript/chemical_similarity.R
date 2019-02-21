setwd("~/GitSoftware/RF-NAPClass/")

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


# read data
data <- read.table("similarity.txt", col.names = c("comp", "pdbid", "tanimoto"))

# remove diagonal were ref is comp
data <- data[!(data$pdbid == data$comp),]
nmr_names <- as.character(read.table("nmr_names.txt")$V1)
xray_names <- as.character(read.table("xray_names.txt")$V1)
data$method <- "other"
data$method[data$pdbid %in% nmr_names] <- "nmr"
data$method[data$pdbid %in% xray_names] <- "xray"
data <- data[(data$pdbid %in% c(nmr_names, xray_names)), ]
rownames(data) <- 1:nrow(data)

# get the maximum tanimoto score
tmp <- plyr::ddply(.data = data, .variables = c("pdbid", "method"), .fun = function(x){data.frame(max_tanimoto = round(max(x$tanimoto), 2), mean_tanimoto = round(mean(x$tanimoto), 2))})
boxplot(max_tanimoto~method, data = tmp, ylim = c(0, 1))


# compare with 
nmr <- plot_clscores(file = "poses_classification_scores_nmr.txt", pred = TRUE, ylab = "CLscore", plotname = "")
xray <- plot_clscores(file = "poses_classification_scores_xray.txt", pred = TRUE, ylab = "CLscore", plotname = "")
all <- rbind(nmr, xray)

tmp <- merge(all, tmp, by = c("pdbid"))
tmp <- tmp[order(tmp$max_tanimoto), ]

tmp_nmr <- tmp[tmp$method=="nmr", c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred", "max_tanimoto", "mean_tanimoto")]
tmp_xray <- tmp[tmp$method=="xray", c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred", "max_tanimoto", "mean_tanimoto")]

write.table(tmp_nmr, file = "similarity_plus_recovery_nmr.txt", quote = FALSE, row.names = FALSE)
write.table(tmp_xray, file = "similarity_plus_recovery_xray.txt", quote = FALSE, row.names = FALSE)


success <- NULL
success_nmr <- NULL
success_xray <- NULL
threshold <- seq(0.3, 1, 0.1)
for (t in threshold){
  which <- tmp$max_tanimoto < t+0.05 & tmp$max_tanimoto > t-0.05
  success <- c(success, mean(tmp$rmsd_pred[which]))
  
  which <- tmp_nmr$max_tanimoto < t+0.05 & tmp_nmr$max_tanimoto > t-0.05
  success_nmr <- c(success_nmr, mean(tmp_nmr$rmsd_pred[which]))
  
  which <- tmp_xray$max_tanimoto < t+0.05 & tmp_xray$max_tanimoto > t-0.05
  success_xray <- c(success_xray, mean(tmp_xray$rmsd_pred[which]))
}
success <- data.frame(threshold, success, success_nmr, success_xray)


plot(success$threshold, success$success, type = "l", lwd = "2", ylim = c(0, 10.0))
lines(success$threshold, success$success_nmr, type = "l", lwd = "2", col = "cyan")
lines(success$threshold, success$success_xray, type = "l", lwd = "2", col = "blue")



