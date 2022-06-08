# GWAS  pipeline#####
#
## Field pea ###
### PSbMV #####
######################
getwd()
# set working directory
# setwd("/group/grains/pulses/adnan/fieldpea/gwas_traitlink")

# R-bundle-Bioconductor/3.10-foss-2019a-R-3.6.1

# load libraries
library(readxl)
library(tidyverse)
library(asreml)
library(lattice)
library(here)
library(janitor)
library(dplyr)
library(stringr)

# dir.create("raw.data")
# dir.create("figures")
# loading data
genofilePath <- "/group/grains/pulses/adnan/fieldpea/f2_geno_filtering/XT_Geno/"
pedPath <- "/group/grains/pulses/adnan/fieldpea/psbmv/"
ozpPath <- "/group/grains/pulses/adnan/fieldpea/fp_psbvm/"

# reading the geno f2 file
geno.data <- readRDS(paste0(genofilePath, "fieldpea_training_candidates_4478IDs_61045SNPs.rds"))
# geno.data[1:5,1:7]
geno.data <- geno.data %>% mutate(Genotype = toupper(Genotype))
geno <- as.matrix(geno.data[, -c(1:2, 4:6)])
# dim(geno)# 8095 61049
# geno[1:5,1:5]
# class(geno.data)

# load pheno file
pbsmv <- read_excel("./raw.data/VirusMETData17_21.xlsx",
  sheet = 3,
  col_types = c("text", "numeric", "numeric", "numeric", "text", "numeric"),
  .name_repair = make_clean_names
)
head(pbsmv)
dim(pbsmv)

# preparing pheno file
pheno <- pbsmv %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(test, as.factor)) %>%
  mutate(across(tray, as.factor)) %>%
  mutate(across(pot, as.factor)) %>%
  rename(line = name) %>%
  mutate(year = case_when(
    r_expt == "PSbMV2017" ~ "2017",
    r_expt == "PSbMV2018" ~ "2018",
    r_expt == "PSbMV2019" ~ "2019", TRUE ~ "2020"
  )) %>%
  select(line, year, everything()) %>%
  group_by(line, year, r_expt) %>%
  summarize(mean_percent_p_sb = mean(percent_p_sb, na.rm = TRUE)) %>%
  mutate(line = toupper(line))

pheno$year <- as.factor(pheno$year)

head(pheno)
dim(pheno) # 476
length(unique(pheno$line)) # 436
table(pheno$r_expt) # identify unique lines in each expt
# PSbMV2017 PSbMV2018 PSbMV2019 PSbMV2020
# 145       102       121       108

list2env(split(pheno, pheno$r_expt), envir = .GlobalEnv) # split data in each expt

# comparing lines between different expt
d.17.18 <- inner_join(PSbMV2017, PSbMV2018, by = "line") # 14 common lines
d.17.19 <- inner_join(PSbMV2017, PSbMV2019, by = "line") # 16
d.17.20 <- inner_join(PSbMV2017, PSbMV2020, by = "line") # 6
d.18.19 <- inner_join(PSbMV2018, PSbMV2019, by = "line") # 3
d.18.20 <- inner_join(PSbMV2018, PSbMV2020, by = "line") # 7
d.19.20 <- inner_join(PSbMV2019, PSbMV2020, by = "line") # 2
# mostly common lines are GREENFEAST, KASPA, PBA WHARTON
# 17 #18 #19 #20
# 17 1
# 18 14  1
# 19 16  3   1
# 20 6   7   2   1

# Extract matching and absent geno lines by comparing with pheno file
com.lines <- geno[geno[, 1] %in% pheno$line, ] # identify the lines genotyped
dim(com.lines) # 424 61046

# prepare a data frame from geno file carrying geno names
lines.sub <- geno[, "Genotype"]
lines.sub.1 <- as.data.frame(lines.sub)
dim(lines.sub.1)
lines.sub.1$lines.sub <- as.character(lines.sub.1$lines.sub)

# subset the lines without the genotypic data (absent from geno file)
abs.lines <- setdiff(pheno$line, lines.sub.1$lines.sub)
abs.lines <- as.data.frame(abs.lines)
dim(abs.lines) # 12

# retain matching lines bw geno and pheno files
genoSUB <- geno[geno[, 1] %in% pheno$line, ]
phenoSUB <- pheno[pheno$line %in% geno[, 1], ]
genoSUB <- genoSUB[order(genoSUB[, 1]), ]
phenoSUB <- phenoSUB[order(phenoSUB$line), ]
# table(phenoSUB$year)

genoSUB <- genoSUB[!duplicated(genoSUB[, 1]), ]
dim(genoSUB) # 424 61045
all(genoSUB[, 1] == unique(phenoSUB$line)) # TRUE
all.equal(genoSUB[, 1], unique(phenoSUB$line)) # TRUE
identical(genoSUB[, 1], unique(phenoSUB$line)) # TRUE
phenoSUB <- phenoSUB[order(phenoSUB$year, phenoSUB$r_expt, phenoSUB$line), ]
dim(phenoSUB) # 462 4


# saving cleaned files
# dir.create("clean_files")
# write.table(phenoSUB, file="./clean_files/PhenoSUB.psbmv.csv", sep = ",", row.names = F, col.names = T)
# dir.create("geno")
# saveRDS(genoSUB, file="./geno/fp_SNPs_linkImputed_psbmv.ATCG.RDS")
# write.table(genoSUB, file="./geno/fp_SNPs_linkImputed_psbmv.ATCG.csv", sep = ",", row.names = F, col.names = T)

# processing geno file for Quality control
PlantID <- genoSUB[, 1]
row.names(genoSUB) <- PlantID
genoSUB <- genoSUB[, -1]
genoSUB <- as.matrix(genoSUB)
genoSUB[1:5, 1:3]
class(genoSUB)
dim(genoSUB) # 424 61045
genoSUB <- apply(genoSUB, 2, as.numeric)
row.names(genoSUB) <- PlantID
length(PlantID) == nrow(genoSUB)

# retain Xs with MAF of higher than 0.05 and lower than 0.95
X <- genoSUB #
maf <- (colSums(genoSUB == 0) + (colSums(genoSUB == 1) * 0.5)) / nrow(genoSUB)
idx <- which(maf >= 0.05 & maf <= 0.95)
X <- genoSUB[, idx]
dim(X) # 428 53536
X[1:3, 1:3] # 424 53507

# preparign the geno file according to GAPIT format
Y <- t(X)
# Y[1:5,1:15]
# options("scipen" = 0, "digits"=9)
Y <- as.data.frame(Y)
Y <- tibble::rownames_to_column(Y, "marker") %>%
  mutate(rs = marker) %>%
  separate(marker, into = c("chrom", "pos")) %>%
  mutate(chrom = factor(chrom,
    levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"),
    labels = c(1, 2, 3, 4, 5, 6, 7)
  )) %>%
  mutate(
    alleles = NA,
    strand = NA,
    assembly = NA,
    center = NA,
    protLSID = NA,
    assayLSID = NA,
    panel = NA,
    QCcode = NA
  ) %>%
  select(rs, alleles, chrom, pos, strand, assembly, center, protLSID, assayLSID, panel, QCcode, everything()) %>%
  mutate(rs = gsub("chr", "rs_chr", rs)) %>%
  rename_at(vars(starts_with("PBA ")), funs(str_replace(., "PBA ", "PBA_")))
# colnames(Y)
# class(Y$pos)
Y$pos <- as.integer(Y$pos)
# Y %>% filter(chrom == "1")

# saveRDS(X, file="./geno/fp_SNPs_linkImputed_psbmv.NUM.MAF.RDS")
# write.table(Y, file="./geno/fp_SNPs_linkImputed_psbmv.NUM.MAF.RDS.hmp.txt", sep = "\t", row.names = F)
# write.table(Y, file="fp_SNPs_linkImputed_psbmv.NUM.MAF.hmp.txt", sep = "\t", row.names = F)

# preparing phenotype file for GAPIT
head(pheno)
pheno.1 <- PSbMV2017 %>%
  mutate(line = gsub("PBA ", "PBA_", line))

list2env(split(phenoSUB, phenoSUB$r_expt), envir = .GlobalEnv) # split data in each expt

pheno.1 <- PSbMV2017 %>%
  full_join(PSbMV2018, by = "line") %>%
  full_join(PSbMV2019, by = "line") %>%
  full_join(PSbMV2020, by = "line") %>%
  select(line, starts_with("mean")) %>%
  rename(
    PSbMV2017 = mean_percent_p_sb.x,
    PSbMV2018 = mean_percent_p_sb.y,
    PSbMV2019 = mean_percent_p_sb.x.x,
    PSbMV2020 = mean_percent_p_sb.y.y
  )

# write.table(pheno.1, file="gwas.pbsmv.pheno.hmp.txt", sep = "\t", row.names = F)

######

# Summary for fieldpea psbmv experiments
fp_psbmv_summary <- data.frame(pheno %>%
  group_by(r_expt) %>%
  summarize(
    count = n(),
    Mean = mean(mean_percent_p_sb, na.rm = TRUE),
    Median = median(mean_percent_p_sb, na.rm = TRUE),
    SD = sd(mean_percent_p_sb, na.rm = TRUE),
    Min. = min(mean_percent_p_sb, na.rm = TRUE),
    Max. = max(mean_percent_p_sb, na.rm = TRUE),
    # CV=sd(mean, na.rm=TRUE)/mean(mean, na.rm=TRUE)*100,
    # St.err= sd(mean, na.rm=TRUE)/sqrt(length(mean))
  )) %>%
  arrange(r_expt)
fp_psbmv_summary <- data.frame(lapply(fp_psbmv_summary, function(y) if (is.numeric(y)) round(y, 2) else y))

fp_psbmv_summary


# First let us visualize the data using boxplots
myboxplot <- function(dataframe, x, y) {
  aaa <- enquo(x)
  bbb <- enquo(y)
  dfname <- enquo(dataframe)
  dataframe %>%
    filter(!is.na(!!aaa), !is.na(!!bbb)) %>%
    # group_by(!! aaa,!! bbb) %>%
    # count() %>%
    ggplot(aes_(fill = aaa, x = aaa, y = bbb)) +
    theme_classic() +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # fill by time point to give different color
    # scale_fill_manual(values = c("", ""))+
    # scale_color_manual(values = c("", ""))
    theme(
      plot.title = element_text(color = "black", size = 12, hjust = 0.5, face = "bold"), # add and modify the title to plot
      axis.title.x = element_text(color = "black", size = 12, face = "bold"), # add and modify title to x axis
      axis.title.y = element_text(color = "black", size = 12, face = "bold")
    ) + # add and modify title to y axis
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    theme(axis.text = element_text(color = "black", size = 10)) + # modify the axis text
    theme(
      legend.title = element_text(colour = "black", size = 16), legend.position = "none",
      legend.text = element_text(colour = "black", size = 14)
    ) + # add and modify the legends
    guides(fill = guide_legend(title = "Environments")) #+
  # stat_summary(fun.y=mean, geom="line",size = 0.1, aes(group=1))  +
  # stat_summary(fun=mean, geom="point",size = 0.1)
}
phenoSUB <- data.frame(phenoSUB)

png("boxplot_all.png", width = 800, height = 600)
p1 <- myboxplot(phenoSUB, x = r_expt, y = mean_percent_p_sb) + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title = "PSbMV 2017-2020", x = "Environments", y = "Disease scale (0-100 %)")
p1
dev.off()

options(bitmapType = "cairo")

# histogram for phenotyped psbmv
png("hist_all.png", width = 800, height = 600)
p2 <- ggplot(data = phenoSUB, aes(x = mean_percent_p_sb)) +
  geom_histogram(binwidth = 10, bins = 10, aes(group = r_expt, fill = r_expt)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~r_expt, ncol = 4) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  xlab("Disease scale (0-100 %)") +
  ylab("Number of lines (n)")
p2
dev.off()


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

###########################
# post gwas anaylsis
############################
# library(devtools)
devtools::install_github("derekmichaelwright/gwaspr", force = TRUE)
library(gwaspr)

# define traits
myTraits <- list_Traits(folder = "GWAS_Results_psbmv/")
myTraits

# define gwas results per trait/exp
myFiles <- list_Result_Files(folder = "GWAS_Results_psbmv/")
myFiles

# getwd()
# setwd( "/group/grains/pulses/adnan/fieldpea/gwas_traitlink/")
myResults_17 <- table_GWAS_Results(
  folder = "GWAS_Results_psbmv/", files = myFiles[[1]],
  threshold = 4, sug.threshold = 3
)

myResults_18 <- table_GWAS_Results(
  folder = "GWAS_Results_psbmv/", files = myFiles[[2]],
  threshold = 4, sug.threshold = 3
)

myResults_19 <- table_GWAS_Results(
  folder = "GWAS_Results_psbmv/", files = myFiles[[3]],
  threshold = 4, sug.threshold = 3
)

myResults_20 <- table_GWAS_Results(
  folder = "GWAS_Results_psbmv/", files = myFiles[[4]],
  threshold = 4, sug.threshold = 3
)

# combining results per trait
myResults <- rbind(myResults_17, myResults_18, myResults_19, myResults_20)
myResults <- myResults %>% arrange(Chromosome, Position)

# selecting markers on chromosome 1 = using previous knowledge
myResults_chr1 <- myResults %>% filter(Chromosome == 1 & `-log10(p)` >= 3.5)

################
# bonferroni threshold
######################
# threshold <- -log10(0.05/nrow(myG))
threshold <- -log10(0.05 / 53508) # 53508 markers
threshold # very stringent so put an arbitrary threshold

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

#################
# list of top markers for each exp
#################
list_Top_Markers_17 <- list_Top_Markers(
  trait = "PSbMV2017", model = "MLM",
  folder = "GWAS_Results/",
  threshold = 4, chroms = 1
)
list_Top_Markers_17$Trait <- "PSbMV2017"

list_Top_Markers_18 <- list_Top_Markers(
  trait = "PSbMV2018", model = "MLM",
  folder = "GWAS_Results/",
  threshold = 4, chroms = 1
)
list_Top_Markers_18$Trait <- "PSbMV2018"

list_Top_Markers_19 <- list_Top_Markers(
  trait = "PSbMV2019", model = "MLM",
  folder = "GWAS_Results/",
  threshold = 4, chroms = 1
)
list_Top_Markers_19$Trait <- "PSbMV2019"

list_Top_Markers_20 <- list_Top_Markers(
  trait = "PSbMV2020", model = "MLM",
  folder = "GWAS_Results/",
  threshold = 3.5, chroms = 1
)
list_Top_Markers_20$Trait <- "PSbMV2020"

list_Top_Markers_all_chr1 <- rbind(list_Top_Markers_17, list_Top_Markers_18, list_Top_Markers_19, list_Top_Markers_20)
list_Top_Markers_all_chr1_unique <- list_Top_Markers_all_chr1[!duplicated(list_Top_Markers_all_chr1[c("POS")]), ]

chr1_unique_snp <- list_Top_Markers_all_chr1_unique$SNP
list_Top_Markers_all_chr1 <- list_Top_Markers_all_chr1 %>%
  left_join(myResults_chr1, by = "SNP") %>%
  select(SNP, Chromosome, Position, Trait.y, "-log10(p).y", effect, pve)

write.csv(list_Top_Markers_all_chr1, "list_Top_Markers_all_chr1.csv", row.names = F)

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



#########################
# plot alleles associated to resistance /susceptible phenotype for psbmv
###################
library(data.table)
options(bitmapType = "cairo")
library(ggplot2)

# subset genotype data for top markers
sub.marker.psbmv.1 <- c()
sub.marker.psbmv <- Y
markers <- as.character(list_Top_Markers_all_chr1$SNP)
for (i in 1:30) {
  sub.marker.psbmv.1[[i]] <- sub.marker.psbmv[sub.marker.psbmv$rs == markers[i], ]
  sub.marker.psbmv.2 <- rbindlist(sub.marker.psbmv.1)
}

# combine with phenotype
sub.marker.psbmv.3 <- data.frame(t(sub.marker.psbmv.2[, -c(2:11)]))
sub.marker.psbmv.3 <- cbind(NAME = row.names(sub.marker.psbmv.3), sub.marker.psbmv.3)
row.names(sub.marker.psbmv.3) <- NULL

sub.marker.psbmv.3 <- sub.marker.psbmv.3 %>%
  full_join(pheno.1, by = c("NAME" = "line")) %>%
  mutate_if(is.factor, as.character)
colnames(sub.marker.psbmv.3)[1:31] <- sub.marker.psbmv.3[1, ]
sub.marker.psbmv.3 <- sub.marker.psbmv.3[-1, ]

# write.csv(sub.marker.psbmv.2. "./geno/geno_psbmv_chr1_marker_only.csv", row.names = F )
# finalise files
sub.marker.psbmv.4 <- sub.marker.psbmv.3[1:424, ]

# to assign factor class to the markers
sub.marker.psbmv.5 <- cbind(id = sub.marker.psbmv.3[1:424, c(1)], data.frame(lapply(sub.marker.psbmv.4[, -1], function(x) as.factor(x))))

# to assign numeric class to the the phenotype data
sub.marker.psbmv.5$PSbMV2017 <- as.numeric(as.character(sub.marker.psbmv.5$PSbMV2017))
sub.marker.psbmv.5$PSbMV2018 <- as.numeric(as.character(sub.marker.psbmv.5$PSbMV2018))
sub.marker.psbmv.5$PSbMV2019 <- as.numeric(as.character(sub.marker.psbmv.5$PSbMV2019))
sub.marker.psbmv.5$PSbMV2020 <- as.numeric(as.character(sub.marker.psbmv.5$PSbMV2020))

# transforming data to long format
sub.marker.psbmv.5_long <- sub.marker.psbmv.5 %>%
  tidyr::gather(key, value, -id, -PSbMV2017, -PSbMV2018, -PSbMV2019, -PSbMV2020)
sub.marker.psbmv.5_long
head(sub.marker.psbmv.5_long)
#                        id PSbMV2017 PSbMV2018 PSbMV2019 PSbMV2020               key value
# 1   05-128-08TGV005#11PS42        40   0.00000        NA        NA rs_chr1.332417863     2
# 2 05H161-06HOS2005-BOG09-2        20   0.00000        NA        NA rs_chr1.332417863     1
# 3   06H246P-08TGV003#11PS5         0        NA        NA        NA rs_chr1.332417863     2
# 4               07H326P004        80  33.33333        NA        NA rs_chr1.332417863     0
# 5            08H016-DO3001        25        NA        NA        NA rs_chr1.332417863     2
# 6            08H201-DO3001       100        NA        NA        NA rs_chr1.332417863     0

# to plot each phenotype individually
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2017, fill = value))+
#  geom_boxplot()+
# facet_wrap(~key)
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2018, fill = value))+
# geom_boxplot()+
#  facet_wrap(~key)
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2019, fill = value))+
# geom_boxplot()+
# facet_wrap(~key)
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2020, fill = value))+
# geom_boxplot()+
# facet_wrap(~key)

# individual plot final
ggplot(sub.marker.psbmv.5_long, aes(x = value, y = PSbMV2017, fill = value)) +
  geom_boxplot() +
  facet_wrap(~key) +
  labs(x = "Alleles", y = "Disease score (0-100 %)", title = "PSbMV2017", fill = "Alleles") +
  geom_rect(
    data = subset(sub.marker.psbmv.5_long, key %in% c(
      "rs_chr1.337689851", "rs_chr1.337756643", "rs_chr1.337765209",
      "rs_chr1.337815258", "rs_chr1.337815367", "rs_chr1.337857927", "rs_chr1.337892153",
      "rs_chr1.337904584", "rs_chr1.337931343", "rs_chr1.337972758", "rs_chr1.338072206",
      "rs_chr1.338076068", "rs_chr1.338109196", "rs_chr1.353954804", "rs_chr1.354002148",
      "rs_chr1.354006151", "rs_chr1.354006590"
    )),
    fill = NA, colour = "red", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

# to plot all plots for phenotype using for loop
exp <- list("PSbMV2017", "PSbMV2018", "PSbMV2019", "PSbMV2020") # create a list of variables to be ploted
for (i in exp) {
  plots <- ggplot(sub.marker.psbmv.5_long, aes(x = value, y = .data[[i]], fill = value)) +
    geom_boxplot() +
    facet_wrap(~key) +
    labs(x = "Alleles", y = "Disease score (0-100 %)", title = i, fill = "Alleles") +
    geom_rect(
      data = subset(sub.marker.psbmv.5_long, key %in% c(
        "rs_chr1.337689851", "rs_chr1.337756643", "rs_chr1.337765209",
        "rs_chr1.337815258", "rs_chr1.337815367", "rs_chr1.337857927", "rs_chr1.337892153",
        "rs_chr1.337904584", "rs_chr1.337931343", "rs_chr1.337972758", "rs_chr1.338072206",
        "rs_chr1.338076068", "rs_chr1.338109196", "rs_chr1.353954804", "rs_chr1.354002148",
        "rs_chr1.354006151", "rs_chr1.354006590"
      )),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  ggsave(plots, filename = paste("alleles res and sus", i, ".pdf"))
}

# GWAS  pipeline##############################################################
#
## Field pea ###
### POWDEY MILDEW #####
######################
getwd()
# set working directory
setwd("/group/grains/pulses/adnan/fieldpea/gwas_traitlink")

# R-bundle-Bioconductor/3.10-foss-2019a-R-3.6.1
# dir.create("raw.data")
dir.create("figures")
# loading data
genofilePath <- "/group/grains/pulses/adnan/fieldpea/f2_geno_filtering/XT_Geno/"
pedPath <- "/group/grains/pulses/adnan/fieldpea/psbmv/"
ozpPath <- "/group/grains/pulses/adnan/fieldpea/fp_psbvm/"

# reading the geno f2 file
geno.data <- readRDS(paste0(genofilePath, "fieldpea_training_candidates_4478IDs_61045SNPs.rds"))
# geno.data[1:5,1:7]
geno.data <- geno.data %>%
  mutate(Genotype = toupper(Genotype)) %>%
  group_by(name = str_replace(name, "PBA ", "PBA_"))
geno <- as.matrix(geno.data[, -c(1:2, 4:6)])
# dim(geno)# 8095 61049
# geno[1:5,1:5]
# class(geno.data)

# lod libraries
library(readxl)
library(tidyverse)
library(asreml)
library(lattice)
library(here)
library(janitor)
library(dplyr)
library(stringr)

# load pheno file
pm <- read_excel("./raw.data/PMHO22 for AgriBio.xlsx",
  sheet = 1,
  col_types = c("text", "numeric", "text", "text", "text"),
  .name_repair = make_clean_names
)
head(pm)
dim(pm) # 300

# the data was recorded as R for resistance, S for susceptible, and Seg for segregating phenotype.
## The data was also differentiated as high/low confidence phenotypes.
# in the working file i have modified
## in the column 5.
## R for low confidence; r for high confidence
## S/s for low confidence; sus for high confidence
## seg for low confidence; suseg for high confidence
## will work only with high confidence data

table(pm$pm_rating)
# R   S   s seg
# 43 138  17 102
table(pm$pm_rating_2)
# R     S     r     s   seg   sus suseg
# 18    25    25    17    86   113    16

# preparing pheno file
pm.1 <- pm %>%
  filter(pm_rating_2 == "r" | pm_rating_2 == "sus" | pm_rating_2 == "suseg") %>%
  mutate(score = recode(pm_rating_2, `r` = "1", `sus` = "0", "suseg" = "2")) %>%
  mutate(name = toupper(name))

head(pm.1)
dim(pm.1) # 476
length(unique(pm.1$name)) # 154
table(pm.1$r_expt) # identify unique lines in each expt
# PPMHO22
# 154
unique(geno.data$Genotype)

# Extract matching and absent geno lines by comparing with pheno file
com.lines.pm <- geno[geno[, 1] %in% pm.1$name, ] # identify the lines genotyped
dim(com.lines.pm) # 148 61046

# prepare a data frame from geno file carrying geno names
lines.sub.pm <- geno[, "Genotype"]
lines.sub.pm <- as.data.frame(lines.sub)
dim(lines.sub.pm)
lines.sub.pm$lines.sub <- as.character(lines.sub.pm$lines.sub)

# subset the lines without the genotypic data (absent from geno file)
abs.lines.pm <- setdiff(pm.1$name, lines.sub.pm$lines.sub)
abs.lines.pm <- as.data.frame(abs.lines)
dim(abs.lines.pm) # 12

# retain matching lines bw geno and pheno files
genoSUB.pm <- geno[geno[, 1] %in% pm.1$name, ]
phenoSUB.pm <- pm.1[pm.1$name %in% geno[, 1], ]
genoSUB.pm <- genoSUB.pm[order(genoSUB.pm[, 1]), ]
phenoSUB.pm <- phenoSUB.pm[order(phenoSUB.pm$name), ]
# table(phenoSUB$year)

genoSUB.pm <- genoSUB.pm[!duplicated(genoSUB.pm[, 1]), ]
dim(genoSUB.pm) # 424 61045
all(genoSUB.pm[, 1] == unique(phenoSUB.pm$name)) # TRUE
all.equal(genoSUB.pm[, 1], unique(phenoSUB.pm$name)) # TRUE
identical(genoSUB.pm[, 1], unique(phenoSUB.pm$name)) # TRUE
# phenoSUB.pm <- phenoSUB.pm[order(phenoSUB.pm$year, phenoSUB.pm$r_expt, phenoSUB.pm$line),]
dim(phenoSUB.pm) # 462 4

phenoSUB.pm <- phenoSUB.pm[, c(3, 6)]

# saving cleaned files
# dir.create("clean_files")
write.table(phenoSUB.pm, file = "./clean_files/PhenoSUB.pm.csv", sep = ",", row.names = F, col.names = T)
write.table(phenoSUB.pm, file = "./clean_files/PhenoSUB.pm.hmp.txt", sep = "\t", row.names = F)

# dir.create("geno")
saveRDS(genoSUB.pm, file = "./geno/fp_SNPs_linkImputed_psbmv.ATCG.RDS")
write.table(genoSUB.pm, file = "./geno/fp_SNPs_linkImputed_psbmv.ATCG.csv", sep = ",", row.names = F, col.names = T)

# processing geno file for Quality control
PlantID <- genoSUB.pm[, 1]
row.names(genoSUB.pm) <- PlantID
genoSUB.pm <- genoSUB.pm[, -1]
genoSUB.pm <- as.matrix(genoSUB.pm)
genoSUB.pm[1:5, 1:3]
class(genoSUB.pm)
dim(genoSUB.pm) # 148 61045
genoSUB.pm <- apply(genoSUB.pm, 2, as.numeric)
row.names(genoSUB.pm) <- PlantID
length(PlantID) == nrow(genoSUB.pm)

# retain Xs with MAF of higher than 0.05 and lower than 0.95
X.pm <- genoSUB.pm #
maf <- (colSums(genoSUB.pm == 0) + (colSums(genoSUB.pm == 1) * 0.5)) / nrow(genoSUB.pm)
idx <- which(maf >= 0.05 & maf <= 0.95)
X.pm <- genoSUB.pm[, idx]
dim(X.pm) # 148 54097
X.pm[1:3, 1:3] # 424 53507

# preparing the geno file according to GAPIT format
Y.pm <- t(X.pm)
# Y[1:5,1:15]
# options("scipen" = 0, "digits"=9)
Y.pm <- as.data.frame(Y.pm)
Y.pm <- tibble::rownames_to_column(Y.pm, "marker") %>%
  mutate(rs = marker) %>%
  separate(marker, into = c("chrom", "pos")) %>%
  mutate(chrom = factor(chrom,
    levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"),
    labels = c(1, 2, 3, 4, 5, 6, 7)
  )) %>%
  mutate(
    alleles = NA,
    strand = NA,
    assembly = NA,
    center = NA,
    protLSID = NA,
    assayLSID = NA,
    panel = NA,
    QCcode = NA
  ) %>%
  select(rs, alleles, chrom, pos, strand, assembly, center, protLSID, assayLSID, panel, QCcode, everything()) %>%
  mutate(rs = gsub("chr", "rs_chr", rs))
# colnames(Y)
# class(Y$pos)
Y.pm$pos <- as.integer(Y.pm$pos)
# Y %>% filter(chrom == "1")

# saveRDS(X, file="./geno/fp_SNPs_linkImputed_psbmv.NUM.MAF.RDS")
# write.table(Y, file="./geno/fp_SNPs_linkImputed_psbmv.NUM.MAF.RDS.hmp.txt", sep = "\t", row.names = F)
write.table(Y.pm, file = "fp_SNPs_linkImputed_pm.NUM.MAF.hmp.txt", sep = "\t", row.names = F)


#########################
# GWAS starting here : Powdery MILDEW - pm
##############################
# use R/4.0.0-foss-2020a in linux/R
# dir.create("gwas_pm")
# setwd("/group/grains/pulses/adnan/fieldpea/gwas_traitlink/gwas_pm/")
# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3", force = TRUE)
library(GAPIT3)

myG <- read.table("fp_SNPs_linkImputed_pm.NUM.MAF.hmp.txt", head = FALSE)
myY <- read.table("PhenoSUB.pm.hmp.txt", head = TRUE)


myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  PCA.total = 3, model = "MLM", SNP.MAF = 0.05
)

###########################
# post gwas anaylsis
############################
# library(devtools)
# devtools::install_github("derekmichaelwright/gwaspr", force = TRUE)
library(gwaspr)
# dir.create("GWAS_Results")
# define traits
myTraits.pm <- list_Traits(folder = "GWAS_Results/")
myTraits.pm

# define gwas results per trait/exp
myFiles.pm <- list_Result_Files(folder = "GWAS_Results/")
myFiles.pm

# getwd()
# setwd( "/group/grains/pulses/adnan/fieldpea/gwas_traitlink/")
myResults_pm <- table_GWAS_Results(
  folder = "GWAS_Results/", files = myFiles[[1]],
  threshold = 4, sug.threshold = 3
)
# combining results per trait
myResults_pm <- myResults_pm %>% arrange(Chromosome, Position)

# selecting markers on chromosome 1 = using previous knowledge
myResults_pm_chr1 <- myResults_pm %>% filter(Chromosome == 1 & `-log10(p)` >= 3.5)

################
# bonferroni threshold
######################
# threshold <- -log10(0.05/nrow(myG))
threshold <- -log10(0.05 / 54098) # 53508 markers
threshold # very stringent so put an arbitrary threshold

########################
# phenotypic variance explained
###########################
# Calc Percent Variance Explained
myResults_pm_chr1$pop.size <- ifelse(grepl("score", myResults_chr1$Trait), "148", "")
myResults_pm_chr1$pop.size <- as.numeric(myResults_chr1$pop.size)

# function pve
getPVE <- function(LOD, N) {
  100 * (1 - 10^((-2 * LOD) / N))
}

myResults_pm_chr1$pve <- getPVE(myResults_pm_chr1$`-log10(p)`, myResults_pm_chr1$pop.size)
# percent variance explained is linked with the LOD values of your peak and can be calculated as such:
# PVE = 1 – 10^((-2*LOD) / n)
# Where n equals the samples size, and LOD is the LOD peak value that can be calculated from the scanone() function.

myResults_pm_chr1 <- myResults_pm_chr1 %>%
  select(SNP, Chromosome, Position, Trait, "-log10(p)", effect, pve)

#################
# list of top markers for each exp
#################
list_Top_Markers_pm <- list_Top_Markers(
  trait = "score", model = "MLM",
  folder = "GWAS_Results/",
  threshold = 3.5, chroms = 1
)
list_Top_Markers_pm$Trait <- "PMHO2022"

library(stringr)
list_Top_Markers_all_pm <- list_Top_Markers_pm %>%
  left_join(myResults_pm_chr1, by = "SNP") %>%
  select(SNP, CHR, POS, Trait.y, "-log10(p).y", effect, pve) %>%
  rename_at(vars(matches(".x|.y")), ~ str_remove(., ".x|.y"))

write.csv(list_Top_Markers_all_pm, "list_Top_Markers_pm.csv", row.names = F)

###################
# Preparing for LDheatmap
####################
myG.ld.pm <- read.table("/group/grains/pulses/adnan/fieldpea/gwas_traitlink/gwas_pm/fp_SNPs_linkImputed_pm.NUM.MAF.hmp.txt", head = F, stringsAsFactors = F)

colnames(myG.ld.pm) <- myG.ld.pm[1, ]
myG.ld.pm <- myG.ld.pm[-1, ]
myG.ld.1.pm <- myG.ld.pm[, c(1, 3, 4, 5, 12:159)]

# rename columns
colnames(myG.ld.1.pm)[1] <- "Name"
colnames(myG.ld.1.pm)[2] <- "Chr"
colnames(myG.ld.1.pm)[3] <- "Pos"
colnames(myG.ld.1.pm)[4] <- "Strand"
myG.ld.1.pm$Strand <- "u"

myG.ld.2.pm <- myG.ld.1.pm[myG.ld.1.pm[, 1] %in% list_Top_Markers_all_pm$SNP, ]
myG.ld.3.pm <- myG.ld.2.pm[, -c(1, 2, 4)]
CEUSNP.pm <- t(myG.ld.3.pm)
colnames(CEUSNP.pm) <- NULL
colnames(CEUSNP.pm) <- CEUSNP.pm[1, ]
CEUSNP.pm <- CEUSNP.pm[-1, ]

snp.dist.pm <- colnames(CEUSNP.pm)
save(snp.dist.pm, file = "snp.dist.pm.RData")
save(CEUSNP.pm, file = "CEUSNP.pm.RData")

################
# LD heatmap POWDERY MILDEW
###############
library(LDheatmap)
require(snpStats)

load("CEUSNP.pm.RData") # load genotype data in numeric format; matrix
load("snp.dist.pm.RData") # marker name as vector

CEUSNP.pm.1 <- apply(CEUSNP.pm, 2, as.numeric) # convert to numeric
gdat_eur.pm <- as(CEUSNP.pm.1, "SnpMatrix")
snp.dist.pm <- as.numeric(snp.dist.pm)
options(bitmapType = "cairo")

# creating plot LD
pdf("rplot.pm.pdf")
myldheatmap.pm <- LDheatmap(gdat_eur.pm, snp.dist.pm, LDmeasure = "r", color = heat.colors(20), SNP.name = c("165050461", "165242176"))
dev.off()

#########################
# plot alleles associated to resistance /susceptible phenotype for psbmv
###################
library(data.table)
options(bitmapType = "cairo")
library(ggplot2)

# subset genotype data for top markers
sub.marker.pm.1 <- c()
sub.marker.pm <- Y.pm
markers.pm <- as.character(list_Top_Markers_all_pm$SNP)
for (i in 1:4) {
  sub.marker.pm.1[[i]] <- sub.marker.pm[sub.marker.pm$rs == markers.pm[i], ]
  sub.marker.pm.2 <- rbindlist(sub.marker.pm.1)
}

# combine with phenotype
sub.marker.pm.3 <- data.frame(t(sub.marker.pm.2[, -c(2:11)]))
sub.marker.pm.3 <- cbind(NAME = row.names(sub.marker.pm.3), sub.marker.pm.3)
row.names(sub.marker.pm.3) <- NULL

# sub.marker.pm.3$NAME <-  as.character(sub.marker.pm.3$NAME)

sub.marker.pm.3 <- sub.marker.pm.3 %>%
  full_join(phenoSUB.pm, by = c("NAME" = "name")) %>%
  mutate(NAME = gsub("PBA ", "PBA_", NAME)) %>%
  filter(!is.na(X1)) %>%
  mutate_if(is.factor, as.character)


colnames(sub.marker.pm.3)[1:5] <- sub.marker.pm.3[1, ]
sub.marker.pm.3 <- sub.marker.pm.3[-1, ]

# write.csv(sub.marker.psbmv.2. "./geno/geno_psbmv_chr1_marker_only.csv", row.names = F )
# finalise files
# sub.marker.pm.4 <- sub.marker.pm.3 %>% filter(!grepl("PBA_", rs))

# to assign factor class to the markers
sub.marker.pm.4 <- cbind(id = sub.marker.pm.3[, c(1)], data.frame(lapply(sub.marker.pm.3[, -1], function(x) as.factor(x))))

# to assign numeric class to the the phenotype data
sub.marker.pm.4$score <- as.numeric(as.character(sub.marker.pm.4$score))

sub.marker.pm.5 <- sub.marker.pm.4
# transforming data to long format
sub.marker.pm.5_long <- sub.marker.pm.5 %>%
  # filter(!is.na(score)) %>%
  tidyr::gather(key, value, -id, -score)
sub.marker.pm.5_long
head(sub.marker.pm.5_long)
#                        id PSbMV2017 PSbMV2018 PSbMV2019 PSbMV2020               key value
# 1   05-128-08TGV005#11PS42        40   0.00000        NA        NA rs_chr1.332417863     2
# 2 05H161-06HOS2005-BOG09-2        20   0.00000        NA        NA rs_chr1.332417863     1
# 3   06H246P-08TGV003#11PS5         0        NA        NA        NA rs_chr1.332417863     2
# 4               07H326P004        80  33.33333        NA        NA rs_chr1.332417863     0
# 5            08H016-DO3001        25        NA        NA        NA rs_chr1.332417863     2
# 6            08H201-DO3001       100        NA        NA        NA rs_chr1.332417863     0

# to plot each phenotype individually
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2017, fill = value))+
#  geom_boxplot()+
# facet_wrap(~key)
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2018, fill = value))+
# geom_boxplot()+
#  facet_wrap(~key)
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2019, fill = value))+
# geom_boxplot()+
# facet_wrap(~key)
# ggplot(sub.marker.psbmv.5_long, aes(x = value , y = PSbMV2020, fill = value))+
# geom_boxplot()+
# facet_wrap(~key)

# individual plot final

ggplot(sub.marker.pm.5_long, aes(x = value)) +
  geom_bar(aes(fill = score)) +
  facet_wrap(~key) +
  labs(x = "Disease rating", y = "Number of lines (n)", title = "PM", fill = "Alleles")
ggsave("alleles res and sus pm.pdf")

#########
# 6. Analysis of variance for haplotypes
## signifcance testing
# Boxplot with mean comparison
# rm(list = ls())
# set working directory
# setwd("C:/Users/user/Desktop/m_download")
# inp.data1 <- read.table("PHHapdata.txt", header = T)
# library(ggplot2)
# library(agricolae)
# library(dplyr)

### hap for combined environment
# value_max = inp.data1 %>% group_by(HAP) %>% summarize(max_value = max(CE))
# hsd=HSD.test(aov(CE~HAP, data=inp.data1), trt = "HAP", group = T)
# hsd
# sig.letters <- hsd$groups[order(row.names(hsd$groups)), ]
# p <- ggplot(inp.data1, aes(x = HAP, y = CE))+
#  geom_boxplot(aes(fill= HAP))+
#  geom_text(data = value_max, aes(x=HAP, y = 0.1 + max_value, label = sig.letters$groups), vjust=0)+
#  stat_boxplot(geom = 'errorbar', width = 0.1)+
#  ggtitle("Haplotype Effect on PH in Combined Environment (CE)") + xlab("") + ylab(""); p
