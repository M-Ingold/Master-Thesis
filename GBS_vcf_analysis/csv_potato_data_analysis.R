library(ape)
library(vcfR)
library(cluster)
library(ggplot2)
library(factoextra)
library(dendextend)
library(readxl)
library(dplyr)
library(ggfortify)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(multcompView)
library(tidyverse)
library(compare)
# 
vcfTetra <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf")
# 
# dosagetetra <- extract.gt(vcfTetra, as.numeric = F)
# 
# dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/0", 0)
# dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/1", 1)
# dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/1/1", 2)
# dosagetetra <- replace(dosagetetra, dosagetetra == "0/1/1/1", 3)
# dosagetetra <- replace(dosagetetra, dosagetetra == "1/1/1/1", 4)
# # to-do: impute NAs
# dosagetetra[is.na(dosagetetra)] <- 0
# #dosagetetra <- as.numeric(na.omit(dosagetetra))
# 
# dosagetetra <- as.numeric(dosagetetra)
# 
# matrix_tetra <- t(matrix(unlist(dosagetetra), nrow = length(dosagetetra)/261, ncol = 261))

# import potato genotype data
GBS02 <- read_excel("../data/GBS2/GBS02_sample-IDs_14.04.21.xlsx", skip = 6)
GBS02 <- rename(GBS02, Identifier = "Sample identifier")
GBS01 <- read_excel("../data/GBS1/22-01-2021_Zuordnung_genotypes_SB.aktualisiert.xlsx")

# poorly merge datasets, maybe do it better later. Too bad!
intersect <- intersect(names(GBS01), names(GBS02))
GBS <- merge(GBS01, GBS02, by=intersect, all = T)

df <- data.frame(Identifier=gsub("Sample_", "",colnames(vcfTetra@gt)[-1]))
Samples <- merge(df, GBS, by="Identifier")

# add sample code to matrix to later substitute proper breed name
#rownames(matrix_tetra) <- Samples$VARIETY

# subset Samples by ADMIXTURE population membership > 90%
admixture=read.table("../scripts/ADMIXTURE/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf.5.Q")
rownames(admixture) <- Samples$VARIETY
admixture <- subset(admixture, V1 > 0.9|V2 > 0.9|V3 > 0.9|V4 > 0.9|V5 > 0.9)

# create dataframe row containing population number
admixture$subpopulation <- NA
for (breed in 1:nrow(admixture)){
  if (admixture[breed,]$V1 > 0.9){
    admixture[breed,]$subpopulation <- "1"
  }
  if (admixture[breed,]$V2 > 0.9){
    admixture[breed,]$subpopulation <- "2"
  }
  if (admixture[breed,]$V3 > 0.9){
    admixture[breed,]$subpopulation <- "3"
  }
  if (admixture[breed,]$V4 > 0.9){
    admixture[breed,]$subpopulation <- "4"
  }
  if (admixture[breed,]$V5 > 0.9){
    admixture[breed,]$subpopulation <- "5"
  }
}

# subset for non-admixed breeds
#matrix_tetra <- matrix_tetra[rownames(admixture), ]



##################################################################################



data <- read_excel("../data/Additional_data/2022-02-14_Screeningdaten_Set_1-7.xlsx", skip = 1)
#sampleNames <- read_excel("../data/Additional_data/22-01-2021_Zuordnung_genotypes_SB.aktualisiert.xlsx")

# Only taking the data for which we have unique identifier of genotyped samples in the 'sampleNames'
data <- data[1:199,]
#sampleNames <- sampleNames[1:189,c(1,5)]

#setdiff(data$`3`, sampleNames$VARIETY)

# Specifying different data frames to concatenate when desired
breedData <- as.data.frame(data[1:12]) # Data common to all measured traits
breedData <- breedData %>% distinct(Variety, .keep_all = T)
breedData <-  replace(breedData, breedData == "DEU", "GER")

merge <- merge(breedData, Samples, by.x = "Variety", by.y = "VARIETY", all=T)

# write all known data of a breed to add parents by hand
write.csv(merge, file="breed_info_full.csv")

# subset varieties with population group and maturity data
#rownames(breedData) <- breedData$Variety
# breedData$Maturity <- as.character(breedData$Maturity)
# x <- merge(breedData, admixture, by = 0)
 


# list of samples with genotyped parents
parentlist <- read_csv("breed_info_full.csv", col_types = "c")

tempparentlist <- parentlist[, c(2,7,8)]
tempparentlist[is.na(tempparentlist)] <- "0"

# if mother or father of a breed is present in the list of genotyped breeds, add it to a list
parents_genotyped <- list()
for (breed in 1:nrow(tempparentlist)){
  if (any(tempparentlist$VARIETY == tempparentlist$Mother[breed])){
    parents_genotyped <- c(parents_genotyped, tempparentlist$VARIETY[breed])
  }
  if (any(tempparentlist$VARIETY == tempparentlist$Father[breed])){
    parents_genotyped <- c(parents_genotyped, tempparentlist$VARIETY[breed])
  }
}

# create a dataframe using the list generated above and add the boolean column to merge dfs
l <- unique(parents_genotyped)
l <- (data.frame(VARIETY = unlist(l)))
l$parent_genotyped <- TRUE

parentlist_new <- merge(parentlist, l, by = "VARIETY", all = T)
parentlist_new$parent_genotyped[is.na(parentlist_new$parent_genotyped)] <- FALSE

write.csv(parentlist_new, file="breed_info_full+genotyped_parents.csv")
