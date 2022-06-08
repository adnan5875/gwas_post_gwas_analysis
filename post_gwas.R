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