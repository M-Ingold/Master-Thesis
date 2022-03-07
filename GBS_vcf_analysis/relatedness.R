library(related)
library(pegas)



# reads into an object of class loci (pegas)
#vcf <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_diploidized.vcf')
vcf <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf', to = 50000)

# convert to genind object (adegenet)
gen <- loci2genind(vcf)
# extract genotype data from genind object as data.frame
dat <- as.data.frame(gen@tab)
# add a column for individuals
# coancestry doesn't like spaces
Samples$VARIETY <- gsub(" ", "_", Samples$VARIETY)

dat <- cbind(indiv=as.character(Samples$VARIETY), dat)
rownames(dat) <- as.character(Samples$VARIETY)
# make sure that column is as.character and not as.factor
dat[[1]] <- as.character(dat[[1]])

relatedness <- coancestry(dat, wang=1, allow.inbreeding = T)

inbreeding <- relatedness$inbreeding

plot(relatedness$relatedness$wang)
plot(inbreeding$LH)

relatedness$relatedness[order(relatedness$relatedness$wang, decreasing = T),]
