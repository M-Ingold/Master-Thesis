library(pegas)
library(ggplot2)
library(reshape2)
library(dplyr)
# required for apparent script
library(outliers)

# parent analysis using apparent script from https://github.com/halelab/apparent
# diplopid VCF in proper format
loci <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf', to = 50000)
#set.seed(1)
# number of SNPs barely effects runtime
#loci <- loci[,sample(ncol(loci), 500)]
analysis_list <- rep("All", 261)
input_df <- as.data.frame(cbind(Samples$VARIETY, analysis_list, loci))
colnames(input_df) <- c("Genotype","key")

#subset <- input_df[1:5, ]
subset <- subset(input_df, Genotype == "ANNABELLE"| Genotype == "Nicola"| Genotype == "MONALISA")

# function from apparent script has to be loaded
apparentOUT <- apparent(subset, MaxIdent=0.10, alpha=0.01, nloci=300, self=TRUE, plot=TRUE, Dyad=TRUE)
