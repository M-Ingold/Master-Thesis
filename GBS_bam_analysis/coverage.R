library(dplyr)

coverage <- read.delim("analysis/coverage")
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
