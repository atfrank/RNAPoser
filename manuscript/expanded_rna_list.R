setwd("~/GitSoftware/RF-NAPClass/")
nmr <- toupper(as.character(read.table("nmr_names.txt")$V1))
xray <- toupper(as.character(read.table("xray_names.txt")$V1))
spa.ln <- toupper(as.character(read.table("spa-ln-training-set.txt")$V1))
rna <- toupper(as.character(read.table("rna.txt")$V1))
spa.ln <- spa.ln[(spa.ln %in% rna)]
spa.ln <- spa.ln[!(spa.ln %in% c(nmr, xray))]

write.table(spa.ln, file = "rna_xray.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
