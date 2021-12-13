library(polyRAD)
library(VariantAnnotation)
library(qqman)
library(pegas)
library(knitr)

# datafull <- VCF2RADdata("data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf", 
#                     possiblePloidies = list(2,4), refgenome = "References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa",
#                     expectedLoci = 35000, expectedAlleles = 35000)


# MAF > 0.05
data005 <- VCF2RADdata("/media/markus/HDD/GBS_Analysis/data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf", 
                        possiblePloidies = list(2,4), refgenome = "/media/markus/HDD/GBS_Analysis/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa",
                        expectedLoci = 7000, expectedAlleles = 7000)


#overdispersionPfull <- TestOverdispersion(datafull, to_test = 1:15)

overdispersionP005 <- TestOverdispersion(data005, to_test = 1:5)

qq(overdispersionP005[["2"]])

qq(overdispersionP005[["3"]]) #best fit

qq(overdispersionP005[["4"]]) 



myhindhe <- HindHe(data005, ploidy = 4)
myhindheByLoc <- colMeans(myhindhe, na.rm = TRUE)

# Hind/He is expected to cluster around (ploidy-1)/ploidy
hist(myhindheByLoc, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus")
abline(v = 0.75, col = "blue", lwd = 2)

# expected HindHe doesn't match the data
ExpectedHindHe(data005, inbreeding = 0, ploidy = 4, overdispersion = 3, reps = 10)

# not used cause simulation and data don't match up
#mean(myhindheByLoc < 0.572)
#keeploci <- names(myhindheByLoc)[myhindheByLoc >= 0.572]
#data <- SubsetByLocus(data, keeploci)

mydataHWE <- IterateHWE(data005, overdispersion = 3)

hist(mydataHWE$alleleFreq, breaks = 20, col = "lightgrey")

mydataPopStruct <- IteratePopStruct(data005, nPcsInit = 15, overdispersion = 3)

myallele <- 1
freqcol <- heat.colors(101)[round(mydataPopStruct$alleleFreqByTaxa[,myallele] * 100) + 1]
plot(mydataPopStruct, pch = 21, bg = freqcol)

plot(mydataPopStruct$ploidyChiSq[1,], mydataPopStruct$ploidyChiSq[2,], 
     xlab = "Chi-squared for diploid model",
     ylab = "Chi-squared for allotetraploid model", log = "xy")
abline(a = 0, b = 1, col = "blue", lwd = 2)

myChiSqRat <- mydataPopStruct$ploidyChiSq[1,] / mydataPopStruct$ploidyChiSq[2,]
myChiSqRat <- tapply(myChiSqRat, mydataPopStruct$alleles2loc, mean)
allelesPerLoc <- as.vector(table(mydataPopStruct$alleles2loc))

library(ggplot2)
ggplot(mapping = aes(x = myhindheByLoc[GetLoci(data005)], y = myChiSqRat, fill = as.factor(allelesPerLoc))) +
  geom_point(shape = 21, size = 3) +
  labs(x = "Hind/He", y = "Ratio of Chi-squared values, diploid to allotetraploid",
       fill = "Alleles per locus") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0.5) +
  scale_fill_brewer(palette = "YlOrRd")

wmgenoPopStruct <- GetWeightedMeanGenotypes(mydataPopStruct)
wmgenoPopStruct[1:10,1:5]



RADdata2VCF(mydataPopStruct, "RADVCF", asSNPs = T, hindhe = F)

Export_Structure(mydataPopStruct, "RADSTRUCTURE")
