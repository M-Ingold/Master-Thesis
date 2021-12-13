library(vcfR)
library(ape)
library(adegenet)
library(poppr)

vcf <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf")

vcf005 <- read.vcfR("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf")
dna <- ape::read.dna("References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta", format = "fasta")
gff <- read.table("References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3", sep="\t", quote="")

chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)
chrom <- proc.chromR(chrom, verbose=T)

chromoqc(chrom)

altAF <- extract.info(vcf,element = "AF", as.numeric = T)
altAF005 <- extract.info(vcf005,element = "AF", as.numeric = T)

png(filename = "allelefrequencies.png", width=1000, height = 500)
par(mfrow=c(1,2))
hist(altAF)
hist(altAF005)
dev.off()

gt <- extract.gt(vcf, element = "GT")

dp <- extract.gt(vcf, element = "DP")

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

