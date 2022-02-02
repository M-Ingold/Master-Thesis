library(dplyr)
library(ggplot2)
#coverage <- read.delim("../data/alignment/samtools_coverage.txt")
#coverage_full <- read.csv("../data/alignment/bedtools_coverage.txt", sep = "\t")
coverage <- read.delim("../data/alignment/coverage_chr01.txt")
coverage <- read.delim("../data/alignment/coverage_chr02.txt")
coverage <- read.delim("../data/alignment/coverage_chr03.txt")


coverage <- subset(coverage, X2 > 14300)

sum(coverage[,3] > 7500)
max(coverage[,3])
min(coverage[,3])

threshold <- quantile(coverage[,3], probs = (0.95))
#hist(coverage$X2, breaks = 100)

mean_depth_density <- ggplot(data = coverage, aes(x = X2)) + 
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) + 
  geom_vline(xintercept = threshold, color = "red") + 
  annotate("text", x = 200000, y = 0.00002, label = paste("F^{-1}~(0.95) == ",threshold), parse = TRUE) + 
  labs(x="Coverage per base", title = "Density plot of coverage per base of Chr01")+
  theme_light()
mean_depth_density


chr01 <- subset(coverage, X.chr == "chr01")
chr02 <- subset(coverage, X.chr == "chr02")

summary(coverage)

plot(coverage[1:10000, ]$pos,coverage[1:10000, ]$coverage)

region = seq(from=12465501,to=12465678)

x <- chr01[chr01$X.chr %in% region, ]

x <- filter(chr01, X.chr %in% region)

x <- coverage %>%
  filter(X.chr == "chr01")

x_region <- x[x$pos > 12465500 | x$pos < 12465679]
