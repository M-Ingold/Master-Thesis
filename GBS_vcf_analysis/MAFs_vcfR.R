# Visualize MAF and PIC

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

# depth per sample per random SNP
test <- extract.gt(vcf, element = "DP", as.numeric = T)[sample(38525,1),]
hist(test,breaks = 50)
#plot(density(test))

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

png(filename = "mafs.png", width=1500, height = 750, res = 150)

par(mfrow=c(1,2))
hist(maf$Frequency, main = "No MAF Filter", xlab = "MAF")
abline(v=mean(maf$Frequency), col="red", lwd=3, lty=2)
text(mean(maf$Frequency)+0.05,length(maf$Frequency)*0.7, labels = round(mean(maf$Frequency),2), col="red")
hist(maf005$Frequency, main = "MAF > 0.05", xlab = "MAF")
abline(v=mean(maf005$Frequency), col="red", lwd=3, lty=2)
text(mean(maf005$Frequency)+0.05,2900, labels = round(mean(maf005$Frequency),2), col="red")

dev.off()



# Polymorphism information content, not used in master thesis
pic <- sapply(maf$Frequency, function(x) 1-(x^2+(1-x)^2)-(2*x^2*(1-x)^2))
pic005 <- sapply(maf005$Frequency, function(x) 1-(x^2+(1-x)^2)-(2*x^2*(1-x)^2))



png(filename = "PIC.png", width=1500, height = 750, res = 150)
par(mfrow=c(1,2))
hist(pic, main = "Histogram of PIC")
abline(v=mean(pic), col="red", lwd=3, lty=2)
text(mean(pic)+0.025,18000, labels = round(mean(pic),2), col="red")
hist(pic005, main = "Histogram of PIC for MAF > 0.05")
abline(v=mean(pic005), col="red", lwd=3, lty=2)
text(mean(pic005)+0.025,1620, labels = round(mean(pic005),2), col="red")
dev.off()

gt <- extract.gt(vcf, element = "GT")

dp <- extract.gt(vcf, element = "DP")
depthoverall <- as.numeric(extract.info(vcf,element = "DP"))
min(depthoverall)
max(depthoverall)
hist(depthoverall)


# testing of some functions

#mydiff <- genetic_diff(vcf) # could have used this for GST

#boolean matrix listing if depth > 61. Use to filter VCF
#bool60 <- matrix(as.numeric(dp) > 61, nrow = 32726, ncol =  264, dimnames = list(rownames(dp), colnames(dp)))

#no polyploidy possible
#my_genind <- vcfR2genind(vcf) 
#gl <- vcfR2genlight(vcf)

#record <- 130
#my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, ref.seq = dna[, 
#                                                                                     gff[record, 4]:gff[record, 5]], start.pos = gff[record, 4], verbose = FALSE)
#my_dnabin1

# convert to genlight: empty object
#x <- vcfR2genlight(vcf005)