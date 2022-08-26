library(dplyr)
library(ggplot2)
library(ggpubr)

# attempt at visually determining a read threshold for VCF filtering later by plotting sequencing depths per base

#coverage_full <- read.csv("../data/alignment/bedtools_coverage.txt", sep = "\t") # too big for my home PC to handle
coverage <- read.delim("../data/alignment/coverage_chr01.txt")
coverage <- read.delim("../data/alignment/coverage_chr02.txt")
coverage <- read.delim("../data/alignment/coverage_chr03.txt")
coverage <- read.delim("../data/alignment/coverage_chr04.txt")
coverage <- read.delim("../data/alignment/coverage_chr05.txt")
coverage <- read.delim("../data/alignment/coverage_chr06.txt")
coverage <- read.delim("../data/alignment/coverage_chr07.txt")
coverage <- read.delim("../data/alignment/coverage_chr08.txt")
coverage <- read.delim("../data/alignment/coverage_chr09.txt")
coverage <- read.delim("../data/alignment/coverage_chr10.txt")
coverage <- read.delim("../data/alignment/coverage_chr11.txt")
coverage <- read.delim("../data/alignment/coverage_chr12.txt")


# filter out read depths that would not pass VCF filters for easier visualization
coverage_sub <- subset(coverage, X4 > 14329)

sum(coverage[,3] > 14329)
max(coverage[,3])
min(coverage[,3])
summary(coverage)


# 95% threshold used in filtering
threshold <- quantile(coverage[,3], probs = (0.95))
threshold_sub <- quantile(coverage_sub[,3], probs = (0.95))

mean_depth_density <- ggplot(data = coverage, aes(x = X4)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = threshold, color = "red") + 
  annotate("text", x = 3000, y = 0.76, label = paste("F^{-1}~(0.95) == ",threshold), parse = TRUE) + 
  labs(x="Coverage per base", title = "All base coverages")+
  theme_bw() +
  xlim(1,5000)
#mean_depth_density


mean_depth_density_sub <- ggplot(data = coverage_sub, aes(x = X4)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = threshold_sub, color = "red") + 
  annotate("text", x = 225000, y = 0.00002, label = paste("F^{-1}~(0.95) == ",threshold_sub), parse = TRUE) + 
  labs(x="Coverage per base", title = "Base coverage > 14329")+
  theme_bw() 
#mean_depth_density_sub

png(filename = "coverage_per_base.png", width = 1500, height = 750, res = 130)
ggarrange(mean_depth_density, mean_depth_density_sub)
dev.off()
