library(vcfR)
library(ape)
library(adegenet)
library(poppr)

vcf <- read.vcfR("data/VCF/freebayes_269_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_tetra.vcf")
dna <- ape::read.dna("References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta", format = "fasta")
gff <- read.table("References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3", sep="\t", quote="")

chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)
chrom <- proc.chromR(chrom, verbose=T)

chromoqc(chrom)

MAF <- extract.info(vcf,element = "AF", as.numeric = T)

hist(MAF)

#no polyploidy possible
#my_genind <- vcfR2genind(vcf) 
#gl <- vcfR2genlight(vcf)

#record <- 130
#my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, ref.seq = dna[, 
#                                                                                     gff[record, 4]:gff[record, 5]], start.pos = gff[record, 4], verbose = FALSE)
#my_dnabin1

my_loci <- vcfR2loci(vcf)