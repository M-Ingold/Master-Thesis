library(vcfR)

x <- read.vcfR("data/VCF/freebayes_269_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf")

info <- vcfR2loci(x)

dim(info)

MAF <- extract.info(x,element = "AF")

hist(MAF)