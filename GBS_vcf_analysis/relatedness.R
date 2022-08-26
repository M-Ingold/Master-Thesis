library(related)
library(pegas)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(AGHmatrix)

#vcfdip <- read.vcf('../data/diploid_VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF_diploidized.vcf', to = 50000)
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
load("Samples")

# coancestry doesn't like spaces
Samples$VARIETY <- gsub(" ", "_", Samples$VARIETY)
# for some reason, sample names are cut off at slashes by coancestry. 
# This means some samples would have the same name and relationships would be overestimated
Samples$VARIETY <- gsub("/", "", Samples$VARIETY)
# uppercase for comparison later
Samples$VARIETY <- toupper(Samples$VARIETY)

dat <- cbind(indiv=as.character(Samples$VARIETY), dat)
rownames(dat) <- as.character(Samples$VARIETY)
# make sure that column is as.character and not as.factor
dat[[1]] <- as.character(dat[[1]])


# pairwise relatedness of samples 
# ritland differs from wang, diploids are the most related there
# trioml takes really long
relatedness <- coancestry(dat, wang=1, 
                          #trioml = 1,
                          #lynchli = 1, lynchrd = 1, ritland = 1, quellergt = 1, dyadml = 1,
                          allow.inbreeding = T, ci95.num.bootstrap = 0)

#save(relatedness, file = "relatedness") # only relevant for saving trioml results
rel <- relatedness$relatedness


parents <- read.csv("../scripts/GBS_vcf_analysis/breed_info_full+genotyped_parents.csv")
parents$VARIETY <- as.character(parents$VARIETY)
parents$Mother <- as.character(parents$Mother)
parents$Father <- as.character(parents$Father)
parents$VARIETY <- gsub(" ", "_", parents$VARIETY)
parents$VARIETY <- gsub("/", "", parents$VARIETY)
parents$Mother <- gsub(" ", "_", parents$Mother)
parents$Mother <- gsub("/", "", parents$Mother)
parents$Father <- gsub(" ", "_", parents$Father)
parents$Father <- gsub("/", "", parents$Father)
parents <- na.omit(parents[, c(2,8,9)])


# VanRaden's method. Comment out if Wang is desired
#####################################
G_FullAutopolyploid <- Gmatrix(matrix_tetra, method="VanRaden", ploidy=4)

# extract inbreeding values from matrix
inbreeding_VR <- data.frame(Nr = 1:length(diag(G_FullAutopolyploid)))
inbreeding_VR$inbreeding <- diag(G_FullAutopolyploid-1) 
inbreeding_VR$sample <- row.names(G_FullAutopolyploid)

# plot inbreeding
png(filename = "VanRaden_inbreeding.png", width = 1200, height = 800, res = 130)

ggplot(inbreeding_VR, aes(Nr, inbreeding, label=sample)) +
  geom_text_repel(max.overlaps = 5, size = 3.5) +
  xlab("") + ylab("Inbreeding coefficient")+
  geom_point(shape=1) +
  theme_bw()

dev.off()

G_FullAutopolyploid_sub <- G_FullAutopolyploid
G_FullAutopolyploid_sub[lower.tri(G_FullAutopolyploid, diag = T)] <- NA

G_matrix_subset <- melt(G_FullAutopolyploid_sub, id.var = rownames(G_FullAutopolyploid_sub)[1])
G_matrix_subset <- na.omit(G_matrix_subset)
ggplot(G_matrix_subset, aes(as.factor(Var1), Var2, group=Var2)) + geom_tile(aes(fill = value)) +
  #geom_text(aes(fill = data$value, label = round(data$value, 1))) +
  scale_fill_gradient(low = "yellow", high = "#D6604D", space = "Lab")

rel <- G_matrix_subset
colnames(rel) <- c("ind1.id","ind2.id","wang")
#####################################

# add column with color, depending on diploids present
rel <- rel %>%
  mutate(color = case_when(
    ind1.id == "INCA_SUN" ~ "red",
    ind1.id == "MAYAN_GOLD" ~ "red",
    ind1.id == "MAYAN_TWILIGHT" ~ "red",
    ind2.id == "INCA_SUN" ~ "red",
    ind2.id == "MAYAN_GOLD" ~ "red",
    ind2.id == "MAYAN_TWILIGHT" ~ "red",
    # too many are marked in blue, not row specific
    # (ind1.id %in% parents$VARIETY & ind2.id %in% parents$Mother) ~ "blue",
    # (ind1.id %in% parents$VARIETY & ind2.id %in% parents$Father) ~ "blue",
    # (ind1.id %in% parents$Mother & ind2.id %in% parents$VARIETY) ~ "blue",
    # (ind1.id %in% parents$Father & ind2.id %in% parents$VARIETY) ~ "blue",
    TRUE ~ "black" # fill rest in black as default
  ))

# parent-offspring pairs
# pretty terrible and slow, hope i can revise this at a later point
for (i in 1:length(rel$ind1.id)) {
  for (j in 1:length(parents$VARIETY)) {
    if (rel$ind1.id[i] == parents$VARIETY[j] & rel$ind2.id[i] == parents$Father[j]) {
      rel$color[i] <- "blue"
    }
    else if (rel$ind1.id[i] == parents$VARIETY[j] & rel$ind2.id[i] == parents$Mother[j]) {
      rel$color[i] <- "blue"
    }
    else if (rel$ind1.id[i] == parents$Mother[j] & rel$ind2.id[i] == parents$VARIETY[j]) {
      rel$color[i] <- "blue"
    }
    else if (rel$ind1.id[i] == parents$Father[j] & rel$ind2.id[i] == parents$VARIETY[j]) {
      rel$color[i] <- "blue"
    }
  }
}
length(rel[rel$color=='blue',]$ind1.id)

# half sibs
# I should really take a tidyverse course
hsib1 <- list()
hsib2 <- list()
n <- 1
for (i in 1:length(parents$VARIETY)) {
  for (j in 1:length(parents$VARIETY)) {
    if (parents$Mother[i] == parents$Mother[j] | parents$Father[i]==parents$Father[j]) {
      if (parents$VARIETY[i] != parents$VARIETY[j]) {
        hsib1[n] <- parents$VARIETY[i]
        hsib2[n] <- parents$VARIETY[j]
        n <- n+1
      }
    }
    else if (parents$Mother[i] == parents$Father[j] | parents$Father[i]==parents$Mother[j]) {
      if (parents$VARIETY[i] != parents$VARIETY[j]) {
        hsib1[n] <- parents$VARIETY[i]
        hsib2[n] <- parents$VARIETY[j]
        n <- n+1
      }
    }
  }
}

for (i in 1:length(rel$ind1.id)) {
  for (j in 1:length(hsib1)) {
    if (rel$ind1.id[i] == hsib1[j] & rel$ind2.id[i] == hsib2[j]) {
      rel$color[i] <- "magenta3"
    }
    else if (rel$ind1.id[i] == hsib2[j] & rel$ind2.id[i] == hsib1[j]) {
      rel$color[i] <- "magenta3"
    }
  }
}
length(rel[rel$color=='magenta3',]$ind1.id)

# full sibling pairs
sib1 <- list()
sib2 <- list()
n <- 1
for (i in 1:length(parents$VARIETY)) {
  for (j in 1:length(parents$VARIETY)) {
    if (parents$Mother[i] == parents$Mother[j] & parents$Father[i]==parents$Father[j]) {
      if (parents$VARIETY[i] != parents$VARIETY[j]) {
        sib1[n] <- parents$VARIETY[i]
        sib2[n] <- parents$VARIETY[j]
        n <- n+1
      }
    }
    else if (parents$Mother[i] == parents$Father[j] & parents$Father[i]==parents$Mother[j]) {
      if (parents$VARIETY[i] != parents$VARIETY[j]) {
        sib1[n] <- parents$VARIETY[i]
        sib2[n] <- parents$VARIETY[j]
        n <- n+1
      }
    }
  }
}

for (i in 1:length(rel$ind1.id)) {
  for (j in 1:length(sib1)) {
    if (rel$ind1.id[i] == sib1[j] & rel$ind2.id[i] == sib2[j]) {
      rel$color[i] <- "green4"
    }
    else if (rel$ind1.id[i] == sib2[j] & rel$ind2.id[i] == sib1[j]) {
      rel$color[i] <- "green4"
    }
  }
}
length(rel[rel$color=='green4',]$ind1.id)



# # for related package results only, not used in master thesis
# inbreeding <- relatedness$inbreeding
# inbreeding$Nr <- 1:261
# #plot(relatedness$relatedness$wang)
# plot(inbreeding$LR)
# inbreeding[order(inbreeding$LH),]
# # sample pairs sorted by Wang's relatedness
# rel_ordered <- relatedness$relatedness[order(relatedness$relatedness$wang, decreasing = T),]
# #relatedness$relatedness[order(relatedness$relatedness$ritland, decreasing = T),]


plot1 <- ggplot(rel, aes(x=rownames(rel), y=wang)) +
  geom_point(shape=1, colour=rel$color) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(y = "relatedness", x='') +
  theme_bw()

related_only <- rel[rel$color=='blue'|rel$color=="green4"|rel$color=="magenta3",]
unrelated <- rel[rel$color=='black'|rel$color=="red",]

plot2 <- ggplot(related_only, aes(x=rownames(related_only), y=wang)) +
  geom_point(shape=1, colour=related_only$color) +
  scale_y_continuous(limits = c(min(rel$wang),max(rel$wang))) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(y = "relatedness", x='') +
  theme_bw()


png(filename = "wang_estimator_of_relatedness.png", width = 1000, height = 500, res = 120)
ggarrange(plot1, plot2, ncol = 2, nrow = 1, labels = "AUTO")
dev.off()



# png(filename = "AGHmatrix_estimator_of_relatedness_MAF_VanRaden.png", width = 1000, height = 500, res = 120)
# ggarrange(plot1, plot2, ncol = 2, nrow = 1, labels = "AUTO")
# dev.off()

# means of relatedness groups
mean(rel[rel$color=='blue'|rel$color=="green4",]$wang)
mean(rel[rel$color=="blue",]$wang)
mean(rel[rel$color=="green4",]$wang)
mean(rel[rel$color=="magenta3",]$wang)
mean(related_only$wang)
mean(unrelated$wang)

# # trioml wasn't used for the master thesis
# png(filename = "trioml_estimator_of_relatedness.png", width = 500, height = 500, res = 150)
# ggplot(relatedness$relatedness, aes(x=pair.no, y=trioml)) +
#   geom_point(shape=1, colour=rel$color)
# dev.off()


# # If all estimators of related are used. Not in master thesis
# df <- melt(rel[,c(1,5:12)] ,  id.vars = c('pair.no','color'))
# 
# png(filename = "seven_estimators_of_relatedness.png", width = 2000, height = 2000, res = 200)
# ggplot(df, aes(pair.no,value)) + 
#   geom_point(shape=1, colour=df$color) + 
#   facet_wrap(~variable, scales = "free")
# dev.off()
# 
# png(filename = "Inbreeding_estimators.png", width = 1200, height = 800, res = 130)
# p1 <- ggplot(inbreeding, aes(Nr, LH, label=ind.id)) + 
#   geom_text_repel() +
#   geom_point(shape=1)
# p2 <- ggplot(inbreeding, aes(Nr, LR, label=ind.id)) + 
#   geom_text_repel() +
#   geom_point(shape=1)
# grid.arrange(p1, p2, nrow = 1)
# dev.off()

