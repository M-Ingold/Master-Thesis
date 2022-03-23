library(related)
library(pegas)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(gridExtra)

# reads into an object of class loci (pegas)
#vcfdip <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf', to = 50000)
# results are roughly the same with and without MAF filtering
#vcf <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf', to = 50000)
vcfMAF <- read.vcf('../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf', to = 50000)


# convert to genind object (adegenet)
#gen <- loci2genind(vcf)
#gendip <- loci2genind(vcfdip)
genMAF <- loci2genind(vcfMAF, ploidy = 4)


# extract genotype data from genind object as data.frame
#dat <- as.data.frame(gendip@tab)
dat <- as.data.frame(genMAF@tab)
#dat <- as.data.frame(gen@tab)

# add a column for individuals
# Sample names from phylogenetic tree script
# coancestry doesn't like spaces
Samples$VARIETY <- gsub(" ", "_", Samples$VARIETY)
# for some reason, sample names are cut off at slashes by coancestry. This means some samples have the same name and relationships are overestimated
Samples$VARIETY <- gsub("/", "", Samples$VARIETY)

dat <- cbind(indiv=as.character(Samples$VARIETY), dat)
rownames(dat) <- as.character(Samples$VARIETY)
# make sure that column is as.character and not as.factor
dat[[1]] <- as.character(dat[[1]])

#dat <- dat[sample(nrow(dat),130),]

# pairwise relatedness of samples 
# ritland differs from wang, diploids are the most related there
# trioml takes too long
relatedness <- coancestry(dat, wang=1, 
                          #trioml = 1, 
                          #lynchli = 1, lynchrd = 1, ritland = 1, quellergt = 1, dyadml = 1,
                          allow.inbreeding = T)

inbreeding <- relatedness$inbreeding
inbreeding$Nr <- 1:261
#plot(relatedness$relatedness$wang)
plot(inbreeding$LR)
inbreeding[order(inbreeding$LH),]
# sample pairs sorted by Wang's relatedness
relatedness$relatedness[order(relatedness$relatedness$wang, decreasing = T),]
#relatedness$relatedness[order(relatedness$relatedness$ritland, decreasing = T),]


png(filename = "wang_estimator_of_relatedness.png", width = 500, height = 500)
ggplot(relatedness$relatedness, aes(x=pair.no, y=wang)) +
  geom_point(shape=1)
dev.off()

# ggplot(relatedness$relatedness, aes(x=pair.no, y=trioml)) +
#   geom_point(shape=1)

df <- melt(relatedness$relatedness[,c(1,6:11)] ,  id.vars = 'pair.no')

png(filename = "Six_estimators_of_relatedness.png", width = 1000, height = 600)
ggplot(df, aes(pair.no,value)) + 
  geom_point(shape=1) + 
  facet_wrap(~variable, scales = "free")
dev.off()

png(filename = "Inbreeding_estimators.png", width = 1000, height = 600)
p1 <- ggplot(inbreeding, aes(Nr, LH, label=ind.id)) + 
  geom_text_repel() +
  geom_point(shape=1)
p2 <- ggplot(inbreeding, aes(Nr, LR, label=ind.id)) + 
  geom_text_repel() +
  geom_point(shape=1)
grid.arrange(p1, p2, nrow = 1)
dev.off()


# clustering using relatedness?
# wang list as matrix
# issue: no comparison with self, probably in different order
# related_matrix <- matrix(relatedness$relatedness$wang, nrow=261, ncol = 261, dimnames = list(Samples$VARIETY,Samples$VARIETY))
# related_wardd2 <- hclust(d= related_matrix, method = "ward.D2")

# calculate relationship matrix
library(AGHmatrix)
G_FullAutopolyploid <- Gmatrix(matrix_tetra, method="Slater", ploidy=4)


