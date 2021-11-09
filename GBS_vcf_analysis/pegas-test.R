library(pegas)

x <- read.vcf('data/VCF/freebayes_269_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf', to = 50000)

X <- as.loci(x)

h <- haplotype(X)