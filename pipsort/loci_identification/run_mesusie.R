#!/usr/bin/env Rscript
  
library(MESuSiE)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
nc <- args[2]
aw <- eval(parse(text=args[3]))
csfile <- args[4]

gwas1 <- read.table("s1_gwas_temp.txt", sep = "\t", header = TRUE, check.names = FALSE, comment.char = "")
gwas2 <- read.table("s2_gwas_temp.txt", sep = "\t", header = TRUE, check.names = FALSE, comment.char = "")

gwas1$Z = gwas1$beta / gwas1$standard_error
gwas2$Z = gwas2$beta / gwas2$standard_error

gwas1 <- gwas1[, c("rsid", "beta", "standard_error", "Z", "n", "pos")]
gwas2 <- gwas2[, c("rsid", "beta", "standard_error", "Z", "n", "pos")]

colnames(gwas1) <- c("SNP", "Beta", "Se", "Z", "N", "POS")
colnames(gwas2) <- c("SNP", "Beta", "Se", "Z", "N", "POS")

rownames(gwas1) <- gwas1$SNP
rownames(gwas2) <- gwas2$SNP

ld1 <- as.matrix(read.table("s1.ld", header = FALSE))
ld2 <- as.matrix(read.table("s2.ld", header = FALSE))

colnames(ld1) = gwas1$SNP
rownames(ld1) = gwas1$SNP
colnames(ld2) = gwas2$SNP
rownames(ld2) = gwas2$SNP

LDlist <- list(eur = ld1, afr = ld2)
gwaslist <- list(eur = gwas1, afr = gwas2)

sink(csfile)
MESuSiE_res<-meSuSie_core(LDlist,gwaslist,L=nc,ancestry_weight=aw)
sink()
sink(outfile)
MESuSiE_res$pip_config
sink()

