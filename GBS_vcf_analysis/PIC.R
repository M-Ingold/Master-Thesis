library(vcfR)


vcf <- read.vcfR("data/VCF/freebayes_269_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf")
maf <- extract.info(vcf,element = "AF", as.numeric = T)

mafgt05 <- maf[maf > 0.05]

pic <- sapply(mafgt05, function(x) 1-(x^2+(1-x)^2)-(2*x^2*(1-x)^2))

png(filename = "PIC_MAF005.png", width=500, height = 500)
hist(pic)
dev.off()