# calculate Tajima's D

library(PopGenome)

# MAF filtered
x <- readData('../data/VCF/test/', format = 'VCF')
x <- neutrality.stats(x)
x@Tajima.D

# not MAF filtered, far more useful
y <- readData('../data/VCF/testfull/', format = 'VCF')
y <- neutrality.stats(y)
y@Tajima.D
