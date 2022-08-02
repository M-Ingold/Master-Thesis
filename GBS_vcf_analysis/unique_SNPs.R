library(ggplot2)

# visualize the number of unique SNPs per sample

snp_counts <- read.delim("../scripts/GBS_vcf_analysis/unique_SNPs_per_sample.txt", skip = 1)

# Rename column for later merge and fix names
names(snp_counts)[names(snp_counts)=="X.3.sample"] <- "Identifier"
snp_counts$Identifier <- gsub("Sample_","",snp_counts$Identifier)

# merge with df containing proper sample names and sort by unique SNP count
load("Samples")
snp_counts_named <- merge(snp_counts, Samples, by = "Identifier")
snp_counts_sorted <- snp_counts_named[order(snp_counts_named$X.6.nHets, decreasing = T),]


#barplot(snp_counts_sorted[snp_counts_sorted$X.6.nHets > 50,]$X.6.nHets, names.arg = snp_counts_sorted[snp_counts_sorted$X.6.nHets > 50,]$VARIETY, las = 2, horiz = T)

png(filename = "unique_SNPs.png", width = 750, height = 500, res = 110)
ggplot(data = snp_counts_sorted[snp_counts_sorted$X.6.nHets > 50,], aes(x = reorder(VARIETY, X.6.nHets), X.6.nHets)) +
  geom_bar(stat="identity") +
  coord_flip() + ylab("Unique SNPs") + xlab("") +
  theme_bw()
dev.off()