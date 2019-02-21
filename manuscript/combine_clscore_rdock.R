setwd("/Users/atfrank/GitSoftware/RF-NAPClass")

get_nslr <- function(tmp, sortby = "score", truth = "status", decreasing = FALSE){
  # get best RMSD and NSLR
  tmp <- tmp[order(tmp[, sortby], decreasing = decreasing), ]
  data.frame(pose=tmp$pose[1], rmsd=tmp$rmsd[1], nslr=nmR::nslr(tmp[, truth]))
}

get_recovery_rates <- function(rmsd, thresholds = seq(1, 3, 0.5)){
  rates <- mean(rmsd)
  # get best RMSD and NSLR
  for (threshold in thresholds){
    rate <- 100*mean(rmsd < threshold)
    rates <- c(rates, rate)
    cat(sprintf("%s %3.1f\n", threshold, rate))
  }
  return(rates)
}

plot_combine <- function(data, systems, ylim = c(-7, 7), which = 1, class_combined = NULL){
  if(is.null(class_combined)){
    # combine rdock and rna-poser scores
    data <- data[(data$id %in% systems), ]
    scales <- seq(0, 200, 1.0)
    rmsds <- NULL
    nslrs <- NULL
    
    # loop over scales
    for (scale in scales){
      data$combined <- -1*scale*log(data$pos_proba+1) + data$score
      combined <- plyr::ddply(.data = data, .variables = c("id"), .fun = get_nslr, sortby = "combined")
      rmsds <- c(rmsds, mean(combined$rmsd))
      nslrs <- c(nslrs, mean(combined$nslr))
    }
    scale_rmsd_nslr <- data.frame(scales, rmsds, nslrs)
    
    # assess using rna-poser classification
    data$combined <- -1*scales[which.min(rmsds)]*log(data$pos_proba+1) + data$score
    combined <- plyr::ddply(.data = data, .variables = c("id"), .fun = get_nslr, sortby = "combined")
    class <- plyr::ddply(.data = data, .variables = c("id"), .fun = get_nslr, sortby = "pos_proba", decreasing = TRUE)
    class_combined <- merge(class, combined, by = c("id"))
    class_combined$rmsd.diff <- class_combined$rmsd.y-class_combined$rmsd.x
    
    # plot 
    class_combined$color <- "forestgreen"
    class_combined$color[class_combined$rmsd.diff > 0] <- "red"
  }
  if (which == 1){
    plot(scales, rmsds, type = "l", lwd = "1.5", xaxt = "n", yaxt = "n", ylab = "RMSD (Å)", ylim = c(2, 6), xlab = "scale (α)")
    labels <- c("0", "50", "100", "150", "200")
    axis(1, at = seq(0, 200, 50), labels = labels, cex.axis = 1.0, col.axis="black", las=2)
    labels <- c("2.0", "3.0", "4.0", "5.0", "6.0")
    axis(2, at = seq(2, 6, 1), labels = labels, cex.axis = 1.0, col.axis="black", las=2)
    abline(v=scales[which.min(rmsds)], lwd ="1.5", lty = "dotted", col ="black")
    abline(h=mean(min(rmsds)), lwd ="1.5", lty = "dotted", col ="black")
    text(75, 2.2, sprintf("min. RMSD: %2.1f Å", mean(min(rmsds))), col = "black")
    text(scales[which.min(rmsds)]-30, 4.5, sprintf("α: %s", scales[which.min(rmsds)]), col = "black")
  } else {
    plot(class_combined$rmsd.diff, type = "h", lwd = "3", ylim = ylim, yaxt = "n", xaxt = "n", col = class_combined$color, pch = 20, xlab = "", ylab = "∆RMSD (Å)", xlim = c(1, nrow(class_combined)))
    #axis(1, at=seq(1, length(systems), 2),labels=systems[seq(1, length(systems), 2)], cex.axis = 1.0, col.axis="black", las=2)
    #axis(3, at=seq(2, length(systems), 2),labels=systems[seq(2, length(systems), 2)], cex.axis = 1.0, col.axis="black", las=2)
    labels <- c("-10.0", "-8.0", "-6.0", "-4.0", "-2.0", "0.0", "2.0", "4.0", "6.0", "8.0", "10.0")
    axis(2, at = seq(-10, 10, 2), labels = labels, cex.axis = 1.0, col.axis="black", las=2)
    #abline(h=0, lty = "dashed")
    rmsd_pos <- mean(class_combined$rmsd.diff[class_combined$rmsd.diff>0])
    rmsd_neg <- mean(class_combined$rmsd.diff[class_combined$rmsd.diff<0])
    abline(h=rmsd_pos, lwd ="1.5", lty = "dotted", col ="red")
    abline(h=rmsd_neg, lwd ="1.5", lty = "dotted", col ="forestgreen")
    #grid(nx = length(systems)-1, ny = 0, col = "lightgray", lty = "dotted")
    text(35, 6, sprintf("∆RMSD>0: %s/%s (%2.1f Å)", sum(class_combined$rmsd.diff>0), nrow(class_combined), abs(rmsd_pos)), col = "red")
    text(35, -8, sprintf("∆RMSD<0: %s/%s (%2.1f Å)", sum(class_combined$rmsd.diff<0), nrow(class_combined), abs(rmsd_neg)), col = "forestgreen")
  }
  return(class_combined)
}

# load data
systems <- read.table("system.txt")$V1
nmr_names <- read.table("nmr_names.txt")$V1
xray_names <- read.table("xray_names.txt")$V1
names <- c("id", "pose", "lig_id", "score", "score.inter", "inter.const", "inter.polar", "inter.rot", "inter.solv", "inter.vdw", "inter.norm", "intra", "dihedral", "intra.dihedral.0", "intra.polar", "intra.polar.0", "intra.repul", "intra.repul.0", "intra.solv", "intra.solv.lig_0", "intra.vdw", "intra.vdw.0", "intra.norm", "restr", "restr.cavity", "restr.norm", "system", "system.const", "system.dihedral", "system.solv", "system.norm", "heavy", "norm", "class", "neg_proba", "pos_proba", "status", "rmsd")
data <- read.table("combined_score_class.txt", col.names = names)

scale <- 1.0
par(fig=(scale)*c(0.0 ,0.5, 0.5, 1.0), new=FALSE)
class_combined <- plot_combine(data, systems, c(-9, 7), which = 1)
par(fig=(scale)*c(0.45, 0.95, 0.5, 1.0), new=TRUE)
class_combined <- plot_combine(data, systems, c(-9, 7), which = 2, class_combined = class_combined)
dev.copy2pdf(file = "combined_score_classifier.pdf")

# get recovery rates for all
get_recovery_rates(class_combined$rmsd.y)
# get recovery rates for NMR
get_recovery_rates(class_combined$rmsd.y[class_combined$id %in% nmr_names])
# get recovery rates for XRAY
get_recovery_rates(class_combined$rmsd.y[class_combined$id %in% xray_names])


# look at interesting cases
class_combined[class_combined$rmsd.diff> 0.7, ]
class_combined[class_combined$rmsd.diff< -1.4, ]


