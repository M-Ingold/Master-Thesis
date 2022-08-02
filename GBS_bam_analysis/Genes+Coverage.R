# Average genes and coverage per window plotted together and tested for correlation
library(ape)
library(zoo)
library(ggplot2)
library(ggpubr)
library(dplyr)

#CDS interpreted as genes
gff <- read.gff("../References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3")
cds <- subset(gff, type == "CDS")
cds <- cds[c(1,4,5)]
cds_sorted <- unique(cds[with(cds, order(seqid, start)), ])
cds_sorted$n <- 1
cds_chr01 <- subset(cds_sorted, seqid == "chr01")

coverage_chr01 <- read.delim("../data/alignment/coverage_chr01.txt", header = F)
coverage_chr02 <- read.delim("../data/alignment/coverage_chr02.txt", header = F)
coverage_chr03 <- read.delim("../data/alignment/coverage_chr03.txt", header = F)
coverage_chr04 <- read.delim("../data/alignment/coverage_chr04.txt", header = F)
coverage_chr05 <- read.delim("../data/alignment/coverage_chr05.txt", header = F)
coverage_chr06 <- read.delim("../data/alignment/coverage_chr06.txt", header = F)
coverage_chr07 <- read.delim("../data/alignment/coverage_chr07.txt", header = F)
coverage_chr08 <- read.delim("../data/alignment/coverage_chr08.txt", header = F)
coverage_chr09 <- read.delim("../data/alignment/coverage_chr09.txt", header = F)
coverage_chr10 <- read.delim("../data/alignment/coverage_chr10.txt", header = F)
coverage_chr11 <- read.delim("../data/alignment/coverage_chr11.txt", header = F)
coverage_chr12 <- read.delim("../data/alignment/coverage_chr12.txt", header = F)
#coverage_all <- read.delim("../data/alignment/bedtools_coverage.txt", header = F)


#sum of genes on a single chromosome
rollingsum1chr <- function(file, window){
  n = 0
  i = 0
  df <- data.frame(start=integer(), end=integer(), sum=integer())
  while (i < max(file$end)) {
    sum <- sum(file[file$start %in% seq(from=i, to=i+window), ]$n)
    df[n, ] <- c(i, i+window, sum)
    n <- n+1
    i <- i+window
  }
  return(df)
}

#sum of genes on all 13 chromosomes
rollingsumallchr <- function(file, window){
  i = 1
  df <- data.frame()
  while(i < 10){
    chr <- c("chr0",i)
    chr <- paste(chr, collapse = "")
    shorteneddf <- subset(file, seqid == chr)
    dftemp <- rollingsum1chr(shorteneddf, window)
    dftemp$chr <- chr
    df <- rbind(df, dftemp)
    i <- i+1
  }
  while(i < 13){
    chr <- c("chr",i)
    chr <- paste(chr, collapse = "")
    shorteneddf <- subset(file, seqid == chr)
    dftemp <- rollingsum1chr(shorteneddf, window)
    dftemp$chr <- chr
    df <- rbind(df, dftemp)
    i <- i+1
  }
  return(df)
}

#chr01genes <- rollingsum1chr(cds_chr01, 1000000)
allgenes <- rollingsumallchr(cds_sorted, 1000000)

#average coverage of the genome 
avgwindow1chr <- function(file, window){
  n = 0
  i = 0
  df <- data.frame(start=integer(), end=integer(), mid=integer(), avg=integer())
  while (i < max(file$V2)) {
    avg <- sum(file[file$V2 %in% seq(from=i, to=i+window), ]$V3)/window
    df[n, ] <- c(i, i+window, (i+i+window)/2, avg)
    n <- n+1
    i <- i+window
    }
  return(df)
}



# I couldn't get the function to call avgwindow1chr properly, so it's calculated beforehand. 
# It causes Error: $ operator is invalid for atomic vectors. WTF?

# reproducible error:
#p <- c("coverage_","chr01")
#file <- paste(p, collapse = "")
#chrcoverage <- avgwindow1chr(file, 1000000)
# works when a big file with all the data is subset and used, file name probably doesn't get recognized as 
# object when used in another function. But: not enough RAM for the whole coverage file

#plotgenerator <- function(chr){
#  p <- c("coverage_",chr)
#  file <- paste(p, collapse = "")
#  chrcoverage <- avgwindow1chr(file, 1000000)
#  sub <- subset(cds_sorted, seqid == chr)
#  plot <- ggplot(sub, aes(start)) +
#            geom_histogram(fill= "blue",alpha=0.5, bins = dim(chrcoverage)[1]) +
#            geom_point(data=chrcoverage, mapping = aes(x=start, y=avg*4)) +
#            geom_smooth(data=chrcoverage, mapping = aes(x=mid, y=avg*4), se=F, color = "black") +
#            scale_y_continuous(name = "Genes per 1 Mb", 
#                               sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Mb")) +
#            ggtitle(chr)
#  print(plot)
#}

#plotgenerator("chr01")



# as the plotgenerator function does not work, the plots had to be generated separately

chr01coveragewindow <- avgwindow1chr(coverage_chr01, 1000000)
chr02coveragewindow <- avgwindow1chr(coverage_chr02, 1000000)
chr03coveragewindow <- avgwindow1chr(coverage_chr03, 1000000)
chr04coveragewindow <- avgwindow1chr(coverage_chr04, 1000000)
chr05coveragewindow <- avgwindow1chr(coverage_chr05, 1000000)
chr06coveragewindow <- avgwindow1chr(coverage_chr06, 1000000)
chr07coveragewindow <- avgwindow1chr(coverage_chr07, 1000000)
chr08coveragewindow <- avgwindow1chr(coverage_chr08, 1000000)
chr09coveragewindow <- avgwindow1chr(coverage_chr09, 1000000)
chr10coveragewindow <- avgwindow1chr(coverage_chr10, 1000000)
chr11coveragewindow <- avgwindow1chr(coverage_chr11, 1000000)
chr12coveragewindow <- avgwindow1chr(coverage_chr12, 1000000)

# correllation tests of genes vs. coverage per bin for all chromosomes 
cor.test(subset(allgenes, chr == "chr01")$sum, chr01coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr02")$sum, chr02coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr03")$sum, chr03coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr04")$sum, chr04coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr05")$sum, chr05coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr06")$sum, chr06coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr07")$sum, chr07coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr08")$sum, chr08coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr08")$sum, chr08coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr09")$sum, chr09coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr10")$sum, chr10coveragewindow$avg[1:60]) # datasets differ by legth for some reason
cor.test(subset(allgenes, chr == "chr11")$sum, chr11coveragewindow$avg)
cor.test(subset(allgenes, chr == "chr12")$sum, chr12coveragewindow$avg)








plot1 <- ggplot(subset(cds_sorted, seqid == "chr01"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 88) +
            geom_point(data=chr01coveragewindow, mapping = aes(x=start, y=avg)) +
            geom_smooth(data=chr01coveragewindow, mapping = aes(x=mid, y=avg), se=F, color = "black") +
            scale_y_continuous(name = "Genes per 1 Mb",
                       sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb")) +
            theme_bw() +
            ggtitle("Chr01") + 
            theme(text = element_text(size = 20))

plot2 <- ggplot(subset(cds_sorted, seqid == "chr02"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 46) +
            geom_point(data=chr02coveragewindow, mapping = aes(x=start, y=avg)) +
            geom_smooth(data=chr02coveragewindow, mapping = aes(x=start, y=avg), se=F, color = "black") +
            theme_bw() +
            scale_y_continuous(name = "Genes per 1 Mb",
                               sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
            ggtitle("Chr02") + 
  theme(text = element_text(size = 20))

plot3 <- ggplot(subset(cds_sorted, seqid == "chr03"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 60) +
            geom_point(data=chr03coveragewindow, mapping = aes(x=start, y=avg)) +
            geom_smooth(data = chr03coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
            theme_bw() +
            scale_y_continuous(name = "Genes per 1 Mb",
                               sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
            ggtitle("Chr03") + 
  theme(text = element_text(size = 20))

plot4 <- ggplot(subset(cds_sorted, seqid == "chr04"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 69) +
  geom_point(data=chr04coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr04coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr04") + 
  theme(text = element_text(size = 20))

plot5 <- ggplot(subset(cds_sorted, seqid == "chr05"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 55) +
  geom_point(data=chr05coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr05coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr05") + 
  theme(text = element_text(size = 20))

plot6 <- ggplot(subset(cds_sorted, seqid == "chr06"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 59) +
  geom_point(data=chr06coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr06coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr06") + 
  theme(text = element_text(size = 20))

plot7 <- ggplot(subset(cds_sorted, seqid == "chr07"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 57) +
  geom_point(data=chr07coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr07coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr07") + 
  theme(text = element_text(size = 20))

plot8 <- ggplot(subset(cds_sorted, seqid == "chr08"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr08coveragewindow)[1]) +
  geom_point(data=chr08coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr08coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr08") + 
  theme(text = element_text(size = 20))

plot9 <- ggplot(subset(cds_sorted, seqid == "chr09"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr09coveragewindow)[1]) +
  geom_point(data=chr09coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr09coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr09") + 
  theme(text = element_text(size = 20))

plot10 <- ggplot(subset(cds_sorted, seqid == "chr10"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr10coveragewindow)[1]) +
  geom_point(data=chr10coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr10coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr10") + 
  theme(text = element_text(size = 20))

plot11 <- ggplot(subset(cds_sorted, seqid == "chr11"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr11coveragewindow)[1]) +
  geom_point(data=chr11coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr11coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr11") + 
  theme(text = element_text(size = 20))

plot12 <- ggplot(subset(cds_sorted, seqid == "chr12"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr12coveragewindow)[1]) +
  geom_point(data=chr12coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr12coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  theme_bw() +
  scale_y_continuous(name = "Genes per 1 Mb",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Mb"))+
  ggtitle("Chr12") + 
  theme(text = element_text(size = 20))


png(filename = 'coveragevsgenedensitychr1-12.png', width = 1500, height = 2400)
ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, ncol = 2, nrow = 6)
dev.off()


# find out the percentage of the genome covered deep enough to find SNPs
# load the subset coverage files and bind them. Not enough RAM to load the full file at once
files <- list.files(path= "../data/alignment",pattern = "coverage_",full.names = T)


DF <- NULL
for (f in files) {
  dat <- read.delim(f, header=F)
  DF <- rbind(DF, dat)
}

coverage <- subset(DF, V3 > 14300)
length(coverage$V1)
hist(coverage$V3)
