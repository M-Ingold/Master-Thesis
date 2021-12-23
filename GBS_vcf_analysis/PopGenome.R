library(PopGenome)

x <- readData('../data/VCF/test/', format = 'VCF')
x <- neutrality.stats(x)

y <- readData('../data/VCF/testfull/', format = 'VCF')
y <- neutrality.stats(y)
