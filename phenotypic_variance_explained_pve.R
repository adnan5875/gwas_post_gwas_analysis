########################
# phenotypic variance explained
###########################
# Calc Percent Variance Explained
myResults_chr1$pop.size <- ifelse(grepl("2017", myResults_chr1$Trait), "143", ifelse(
  grepl("2018", myResults_chr1$Trait), "101", ifelse(
    grepl("_2019", myResults_chr1$Trait), "115", "103"
  )
))
myResults_chr1$pop.size <- as.numeric(myResults_chr1$pop.size)

# function pve
getPVE <- function(LOD, N) {
  100 * (1 - 10^((-2 * LOD) / N))
}

myResults_chr1$pve <- getPVE(myResults_chr1$`-log10(p)`, myResults_chr1$pop.size)
# percent variance explained is linked with the LOD values of your peak and can be calculated as such:
# PVE = 1 – 10^((-2*LOD) / n)
# Where n equals the samples size, and LOD is the LOD peak value that can be calculated from the scanone() function.

myResults_chr1 <- myResults_chr1 %>%
  select(SNP, Chromosome, Position, Trait, "-log10(p)", effect, pve)
