#########################
# GWAS starting here
##############################
# use R/4.0.0-foss-2020a in linux/R
# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3", force = TRUE)
library(GAPIT3)

myG <- read.table("fp_SNPs_linkImputed_psbmv.NUM.MAF.hmp.txt", head = FALSE)
myY <- read.table("gwas.psbmv.pheno.hmp.txt", head = TRUE)

myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  PCA.total = 3, model = "MLM", SNP.MAF = 0.1
)
