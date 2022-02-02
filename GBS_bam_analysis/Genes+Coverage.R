#Average genes and coverage per window tested for correlation
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
#coverage_all <- read.delim("../data/alignment/samtools_coverage.txt", header = F)


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
  #df <- data.frame(start=integer(), end=integer(), sum=integer(), chr=character())
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
#allgenes <- rollingsumallchr(cds_sorted, 1000000)

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

#different results, why?
# avgwindow1chrtest <- function(file, window){
#   df <- data.frame(matrix(NA, nrow = length(seq(0, max(file$V2), by=window))-1, ncol = 2))
#   
#   df$avg <- tapply(file$V3, cut(file$V2, seq(0, max(file$V2), by=window)), mean)
#   
#   df$start <- seq(0, max(file$V2), by=window)[1:length(seq(0, max(file$V2), by=window))-1]
#   
#   return(df)
# }
# 
# avgwindow1chrtest(coverage_chr01, 1000000)

# chr01coverage <- avgwindow1chr(coverage_chr01, 1000000)
# 
# plot(chr01genes$sum, chr01coverage$avg)
# 
# abline(lm(chr01coverage$avg ~ chr01genes$sum))
# 
# cor.test(chr01genes$sum, chr01coverage$avg)
# 
# 
# chr02coverage <- avgwindow1chr(coverage_chr02, 1000000)
# plot(subset(allgenes, chr == "chr02")$sum, chr02coverage[1:460,]$avg)
# 
# abline(lm(chr02coverage[1:460,]$avg ~ subset(allgenes, chr == "chr02")$sum))
# 
# cor.test(subset(allgenes, chr == "chr02")$sum, chr02coverage[1:460,]$avg)
# 
# 
# 
# chr03coverage <- avgwindow1chr(coverage_chr03, 1000000)
# plot(subset(allgenes, chr == "chr03")$sum, chr03coverage$avg)
# 
# abline(lm(chr03coverage$avg ~ subset(allgenes, chr == "chr03")$sum))
# 
# cor.test(subset(allgenes, chr == "chr03")$sum, chr03coverage$avg)

# #display gene density
# d <- hist(subset(cds_sorted, seqid == "chr01")$start, 
#           xlim = c(0,max(subset(cds_sorted, seqid == "chr01")$start)),
#           breaks = 40)
# #density plot
# d <- plot(chr01coverage$avg, add = T)
# d

# rollwind <- rollmean(coverage_chr01$V3, 1000000)
# plot(rollwind)
# plot(chr01genes$sum)
# plot(chr02coverage$avg)

# plot1 <- ggplot(subset(cds_sorted, seqid == "chr01"), aes(start)) +
#             geom_histogram(fill= "blue",alpha=0.5, bins = 88) +
#             geom_point(data=chr01coverage, mapping = aes(x=start, y=avg*4)) +
#             geom_smooth(data=chr01coverage, mapping = aes(x=mid, y=avg*4), se=F, color = "black") +
#             scale_y_continuous(name = "Genes per bin", 
#                        sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Gb")) +
#             ggtitle("Chr01")
# 
# plot2 <- ggplot(subset(cds_sorted, seqid == "chr02"), aes(start)) +
#             geom_histogram(fill= "blue",alpha=0.5, bins = 46) +
#             geom_point(data=chr02coverage, mapping = aes(x=start, y=avg*4)) +
#             geom_smooth(data=chr02coverage, mapping = aes(x=start, y=avg*4), se=F, color = "black") +
#             #theme_bw() +
#             scale_y_continuous(name = "Genes per bin", 
#                                sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Gb"))+
#             ggtitle("Chr02")
# 
# plot3 <- ggplot(subset(cds_sorted, seqid == "chr03"), aes(start)) +
#             geom_histogram(fill= "blue",alpha=0.5, bins = 46) +
#             geom_point(data=chr03coverage, mapping = aes(x=start, y=avg*4)) +
#             geom_smooth(data = chr03coverage, mapping = aes(x=start, y=avg*4), se = F, color = "black")+
#             #geom_smooth(data=chr02coverage, mapping = aes(x=start, y=avg*4), se=T, color = "black") +
#             #theme_ipsum() +
#             scale_y_continuous(name = "Genes per bin", 
#                                sec.axis = sec_axis(trans=~./5, name = "Average coverage per 1Gb"))+
#             ggtitle("Chr03")
# 
# theme_update(text = element_text(size=20))
# png(filename = 'coveragevsgenedensitychr1-3.png', width = 1000, height = 1000)
# ggarrange(plot1, plot2, plot3)
# dev.off()


# I couldn't get the function to call avgwindow1chr properly, so it's calculated beforehand. 
# It causes Error: $ operator is invalid for atomic vectors. WTF?

# p <- c("coverage_","chr01")
# file <- paste(p, collapse = "")
# chrcoverage <- avgwindow1chr(file, 1000000)

# works when a big file with all the data is subset and used, file name doesn't get recognized as object when used in another function
plotgenerator <- function(chr){
  p <- c("coverage_",chr)
  file <- paste(p, collapse = "")
  chrcoverage <- avgwindow1chr(file, 1000000)
  sub <- subset(cds_sorted, seqid == chr)
  plot <- ggplot(sub, aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = dim(chrcoverage)[1]) +
            geom_point(data=chrcoverage, mapping = aes(x=start, y=avg*4)) +
            geom_smooth(data=chrcoverage, mapping = aes(x=mid, y=avg*4), se=F, color = "black") +
            scale_y_continuous(name = "Genes per bin", 
                               sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Gb")) +
            ggtitle(chr)
  print(plot)
}

plotgeneratortest <- function(chr){
  p <- c("coverage_",chr)
  file <- paste(p, collapse = "")
  #print(typeof(file))
  chrcoverage <- avgwindow1chr(subset(coverage_all, V1 == "chr01"), 1e6)
  sub <- subset(cds_sorted, seqid == chr)
  plot <- ggplot(sub, aes(start)) +
    geom_histogram(fill= "blue",alpha=0.5, bins = dim(chrcoverage)[1]) +
    geom_point(data=chrcoverage, mapping = aes(x=start, y=avg*4)) +
    geom_smooth(data=chrcoverage, mapping = aes(x=start, y=avg*4), se=F, color = "black") +
    scale_y_continuous(name = "Genes per bin", 
                       sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Gb")) +
    ggtitle(chr)
  print(plot)
}

plotgeneratortest("chr01")
# 
# plot01 <- plotgenerator(chr = "chr01")
# plot02 <- plotgenerator("chr02")
# plot03 <- plotgenerator("chr03")
# plot04 <- plotgenerator("chr04")
# plot05 <- plotgenerator("chr05")
# plot06 <- plotgenerator("chr06")
# plot07 <- plotgenerator("chr07")
# plot08 <- plotgenerator("chr08")
# plot09 <- plotgenerator("chr09")
# plot10 <- plotgenerator("chr10")
# plot11 <- plotgenerator("chr11")
# plot12 <- plotgenerator("chr12")

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


plot1 <- ggplot(subset(cds_sorted, seqid == "chr01"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 88) +
            geom_point(data=chr01coveragewindow, mapping = aes(x=start, y=avg)) +
            geom_smooth(data=chr01coveragewindow, mapping = aes(x=mid, y=avg), se=F, color = "black") +
            scale_y_continuous(name = "Genes per bin",
                       sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb")) +
            ggtitle("Chr01")

plot2 <- ggplot(subset(cds_sorted, seqid == "chr02"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 46) +
            geom_point(data=chr02coveragewindow, mapping = aes(x=start, y=avg)) +
            geom_smooth(data=chr02coveragewindow, mapping = aes(x=start, y=avg), se=F, color = "black") +
            #theme_bw() +
            scale_y_continuous(name = "Genes per bin",
                               sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
            ggtitle("Chr02")

plot3 <- ggplot(subset(cds_sorted, seqid == "chr03"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 60) +
            geom_point(data=chr03coveragewindow, mapping = aes(x=start, y=avg)) +
            geom_smooth(data = chr03coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
            #theme_ipsum() +
            scale_y_continuous(name = "Genes per bin",
                               sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
            ggtitle("Chr03")

plot4 <- ggplot(subset(cds_sorted, seqid == "chr04"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 69) +
  geom_point(data=chr04coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr04coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr04")

plot5 <- ggplot(subset(cds_sorted, seqid == "chr05"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 55) +
  geom_point(data=chr05coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr05coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr05")

plot6 <- ggplot(subset(cds_sorted, seqid == "chr06"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 59) +
  geom_point(data=chr06coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr06coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr06")

plot7 <- ggplot(subset(cds_sorted, seqid == "chr07"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = 57) +
  geom_point(data=chr07coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr07coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr07")

plot8 <- ggplot(subset(cds_sorted, seqid == "chr08"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr08coveragewindow)[1]) +
  geom_point(data=chr08coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr08coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr08")

plot9 <- ggplot(subset(cds_sorted, seqid == "chr09"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr09coveragewindow)[1]) +
  geom_point(data=chr09coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr09coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr09")

plot10 <- ggplot(subset(cds_sorted, seqid == "chr10"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr10coveragewindow)[1]) +
  geom_point(data=chr10coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr10coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr10")

plot11 <- ggplot(subset(cds_sorted, seqid == "chr11"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr11coveragewindow)[1]) +
  geom_point(data=chr11coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr11coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr11")

plot12 <- ggplot(subset(cds_sorted, seqid == "chr12"), aes(start)) +
  geom_histogram(fill= "blue",alpha=0.5, bins = dim(chr12coveragewindow)[1]) +
  geom_point(data=chr12coveragewindow, mapping = aes(x=start, y=avg)) +
  geom_smooth(data = chr12coveragewindow, mapping = aes(x=start, y=avg), se = F, color = "black")+
  #theme_ipsum() +
  scale_y_continuous(name = "Genes per bin",
                     sec.axis = sec_axis(trans=~./1, name = "Average coverage per 1Gb"))+
  ggtitle("Chr12")


theme_update(text = element_text(size=20))
png(filename = 'coveragevsgenedensitychr1-12.png', width = 1500, height = 2500)
ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, ncol = 2, nrow = 6)
dev.off()


