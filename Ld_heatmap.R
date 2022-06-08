###################
# Preparing for LDheatmap
####################

myG.ld <- read.table("fp_SNPs_linkImputed_psbmv.NUM.MAF.hmp.txt", head = F, stringsAsFactors = F)
class(myG.ld)
dim(myG.ld)
myG.ld[1:5, 1:17]
colnames(myG.ld) <- myG.ld[1, ]
myG.ld <- myG.ld[-1, ]

myG.ld.1 <- myG.ld[, c(1, 3, 4, 5, 12:435)]
myG.ld.1[1:5, 1:17]

# rename columns
colnames(myG.ld.1)[1] <- "Name"
colnames(myG.ld.1)[2] <- "Chr"
colnames(myG.ld.1)[3] <- "Pos"
colnames(myG.ld.1)[4] <- "Strand"
myG.ld.1$Strand <- "u"

myG.ld.2 <- myG.ld.1[myG.ld.1[, 1] %in% chr1_unique_snp, ]
dim(myG.ld.3)
myG.ld.2[1:5, 1:12]
myG.ld.3 <- myG.ld.2[, -c(1, 2, 4)]
myG.ld.3[1:5, 1:12]
dim()

CEUSNP <- t(myG.ld.2[, c(3, 5:428)])
CEUSNP[1:17, 1:17]
colnames(CEUSNP) <- NULL
colnames(CEUSNP) <- CEUSNP[1, ]
CEUSNP <- CEUSNP[-1, ]

snp.dist <- colnames(CEUSNP)
save(snp.dist, file = "snp.dist.RData")
save(CEUSNP, file = "CEUSNP.RData")

################
# LD heatmap
###############
library(LDheatmap)
require(snpStats)

load("CEUSNP.RData") # load genotype data in numeric format; matrix
load("snp.dist.RData") # marker name as vector

CEUSNP.1 <- apply(CEUSNP, 2, as.numeric) # convert to numeric
gdat_eur <- as(CEUSNP.1, "SnpMatrix")

snp.dist <- as.numeric(snp.dist)

# creating plot LD
pdf("rplot.pdf")
myldheatmap <- LDheatmap(gdat_eur, snp.dist,
  LDmeasure = "r", color = heat.colors(20),
  SNP.name = c("332417863", "354006590"), text = FALSE
)
dev.off()
