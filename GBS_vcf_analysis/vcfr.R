library(vcfR)
library(ape)
library(adegenet)
library(poppr)

vcf <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf")

vcf005 <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf")
# dna <- ape::read.dna("References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta", format = "fasta")
# gff <- read.table("References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3", sep="\t", quote="")
# 
# chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)
# chrom <- proc.chromR(chrom, verbose=T)
# 
# chromoqc(chrom)

altAF <- extract.info(vcf,element = "AF", as.numeric = T)
altAF005 <- extract.info(vcf005,element = "AF", as.numeric = T)

png(filename = "allelefrequencies.png", width=1000, height = 500)
par(mfrow=c(1,2))
hist(altAF, main = "Histogram of alternative allele frequency")
hist(altAF005, main = "Histogram of alternative allele frequency for MAF > 0.05")
dev.off()

maf <- maf(vcf,2)
maf <- as.data.frame(maf)

maf005 <- maf(vcf005,2)
maf005 <- as.data.frame(maf005)

png(filename = "mafs.png", width=1000, height = 500)
par(mfrow=c(1,2))
hist(maf$Frequency, main = "Histogram of MAF")
hist(maf005$Frequency, main = "Histogram of MAF > 0.05")
dev.off()


hist(test$Frequency)
sum(test$Frequency >= 0.05)
sum(test$Frequency >= 0.05)/nrow(test)

# PIC
pic <- sapply(maf$Frequency, function(x) 1-(x^2+(1-x)^2)-(2*x^2*(1-x)^2))
pic005 <- sapply(maf005$Frequency, function(x) 1-(x^2+(1-x)^2)-(2*x^2*(1-x)^2))



png(filename = "PIC.png", width=1000, height = 500)
par(mfrow=c(1,2))
hist(pic, main = "Histogram of PIC")
abline(v=mean(pic), col="red", lwd=3, lty=2)
hist(pic005, main = "Histogram of PIC for MAF > 0.05")
abline(v=mean(pic005), col="red", lwd=3, lty=2)
dev.off()

gt <- extract.gt(vcf, element = "GT")

dp <- extract.gt(vcf, element = "DP")
depthoverall <- as.numeric(extract.info(vcf,element = "DP"))
min(depthoverall)
max(depthoverall)
hist(depthoverall)

#population list needed
#mydiff <- genetic_diff(vcf)

#boolean matrix listing if depth > 61. Use to filter VCF
bool60 <- matrix(as.numeric(dp) > 61, nrow = 32726, ncol =  264, dimnames = list(rownames(dp), colnames(dp)))

#no polyploidy possible
#my_genind <- vcfR2genind(vcf) 
#gl <- vcfR2genlight(vcf)

#record <- 130
#my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, ref.seq = dna[, 
#                                                                                     gff[record, 4]:gff[record, 5]], start.pos = gff[record, 4], verbose = FALSE)
#my_dnabin1

# convert to genlight: empty object
#x <- vcfR2genlight(vcf005)