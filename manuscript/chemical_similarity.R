setwd("~/Documents/GitHub/RNAPoser/manuscript/")

# read data
data <- read.table("chemical_similarity.txt", col.names = c("pdbid", "comp", "tanimoto"))
# remove diagonal were ref is comp
data$pdbid <- as.character(data$pdbid)
data$comp <- as.character(data$comp)
data <- data[!(data$pdbid == data$comp),]
# get summary statistics
loo <- read.table(file = "poses_classification_scores_xray_loo.txt", col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))
valid <- read.table(file = "poses_classification_scores_xray_valid.txt", col.names = c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred"))
all <- rbind(loo, valid)

loo_names <- toupper(as.character(loo$pdbid))
valid_names <- toupper(as.character(valid$pdbid))

data$method <- "other"
data$method[data$pdbid %in% loo_names] <- "loo"
data$method[data$pdbid %in% valid_names] <- "valid"
data <- data[(data$pdbid %in% c(loo_names, valid_names)), ]
rownames(data) <- 1:nrow(data)

# get the maximum tanimoto score
tmp <- plyr::ddply(.data = data, .variables = c("pdbid", "method"), .fun = function(x){data.frame(max_tanimoto = round(max(x$tanimoto), 2), mean_tanimoto = round(mean(x$tanimoto), 2))})
boxplot(max_tanimoto~method, data = tmp, ylim = c(0, 1))


# compare with 
tmp <- merge(all, tmp, by = c("pdbid"))
tmp <- tmp[order(tmp$max_tanimoto), ]

tmp_loo <- tmp[tmp$method=="loo", c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred", "max_tanimoto", "mean_tanimoto")]
tmp_valid <- tmp[tmp$method=="valid", c("pdbid", "cs_true", "rmsd_true", "cs_pred", "rmsd_pred", "max_tanimoto", "mean_tanimoto")]

write.table(tmp_loo, file = "similarity_plus_recovery_loo.txt", quote = FALSE, row.names = FALSE)
write.table(tmp_valid, file = "similarity_plus_recovery_valid.txt", quote = FALSE, row.names = FALSE)

cor(tmp_valid$rmsd_pred, tmp_valid$max_tanimoto, method = "kendall")
cor(tmp_valid$rmsd_pred, tmp_valid$mean_tanimoto, method = "kendall")


success <- NULL
success_loo <- NULL
success_valid <- NULL
incr <- 0.05
incr2 <- incr/2
threshold <- seq(0.3, 1, incr)
for (t in threshold){
  lower <- t+incr2
  upper <- t-incr2
  which <- tmp$max_tanimoto < lower & tmp$max_tanimoto > upper
  success <- c(success, mean(tmp$rmsd_pred[which], na.rm = TRUE))
  
  which <- tmp_loo$max_tanimoto < lower & tmp_loo$max_tanimoto > upper
  success_loo <- c(success_loo, mean(tmp_loo$rmsd_pred[which], na.rm = TRUE))
  
  which <- tmp_valid$max_tanimoto < lower & tmp_valid$max_tanimoto > upper
  success_valid <- c(success_valid, mean(tmp_valid$rmsd_pred[which], na.rm = TRUE))
}
success <- data.frame(threshold, success, success_loo, success_valid)


plot(success$threshold, success$success, type = "s", lwd = "2", ylim = c(0, 10.0))
lines(success$threshold, success$success_loo, type = "s", lwd = "2", col = "cyan")
lines(success$threshold, success$success_valid, type = "s", lwd = "2", col = "blue")



