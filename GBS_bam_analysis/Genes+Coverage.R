#Average genes and coverage per window tested for correlation
library(ape)
library(zoo)
library(ggplot2)
library(ggpubr)
#CDS interpreted as genes
gff <- read.gff("References/DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3")
cds <- subset(gff, type == "CDS")
cds <- cds[c(1,4,5)]
cds_sorted <- unique(cds[with(cds, order(seqid, start)), ])
cds_sorted$n <- 1
cds_chr01 <- subset(cds_sorted, seqid == "chr01")

coverage_chr01 <- read.delim("data/alignment/coverage_chr01.txt", header = F)
coverage_chr02 <- read.delim("data/alignment/coverage_chr02.txt", header = F)
coverage_chr03 <- read.delim("data/alignment/coverage_chr03.txt", header = F)


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

chr01genes <- rollingsum1chr(cds_chr01, 1000000)
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

chr01coverage <- avgwindow1chr(coverage_chr01, 1000000)

plot(chr01genes$sum, chr01coverage$avg)

abline(lm(chr01coverage$avg ~ chr01genes$sum))

cor.test(chr01genes$sum, chr01coverage$avg)


chr02coverage <- avgwindow1chr(coverage_chr02, 1000000)
plot(subset(allgenes, chr == "chr02")$sum, chr02coverage[1:460,]$avg)

abline(lm(chr02coverage[1:460,]$avg ~ subset(allgenes, chr == "chr02")$sum))

cor.test(subset(allgenes, chr == "chr02")$sum, chr02coverage[1:460,]$avg)



chr03coverage <- avgwindow1chr(coverage_chr03, 1000000)
plot(subset(allgenes, chr == "chr03")$sum, chr03coverage$avg)

abline(lm(chr03coverage$avg ~ subset(allgenes, chr == "chr03")$sum))

cor.test(subset(allgenes, chr == "chr03")$sum, chr03coverage$avg)

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

plot1 <- ggplot(subset(cds_sorted, seqid == "chr01"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 88) +
            geom_point(data=chr01coverage, mapping = aes(x=start, y=avg*4)) +
            geom_smooth(data=chr01coverage, mapping = aes(x=mid, y=avg*4), se=F, color = "black") +
            scale_y_continuous(name = "Genes per bin", 
                       sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Gb")) +
            ggtitle("Chr01")

plot2 <- ggplot(subset(cds_sorted, seqid == "chr02"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 46) +
            geom_point(data=chr02coverage, mapping = aes(x=start, y=avg*4)) +
            geom_smooth(data=chr02coverage, mapping = aes(x=start, y=avg*4), se=F, color = "black") +
            #theme_bw() +
            scale_y_continuous(name = "Genes per bin", 
                               sec.axis = sec_axis(trans=~./4, name = "Average coverage per 1Gb"))+
            ggtitle("Chr02")

plot3 <- ggplot(subset(cds_sorted, seqid == "chr03"), aes(start)) +
            geom_histogram(fill= "blue",alpha=0.5, bins = 46) +
            geom_point(data=chr03coverage, mapping = aes(x=start, y=avg*5)) +
            geom_smooth(data = chr03coverage, mapping = aes(x=start, y=avg*5), se = F, color = "black")+
            #geom_smooth(data=chr02coverage, mapping = aes(x=start, y=avg*8), se=T, color = "black") +
            #theme_ipsum() +
            scale_y_continuous(name = "Genes per bin", 
                               sec.axis = sec_axis(trans=~./5, name = "Average coverage per 1Gb"))+
            ggtitle("Chr03")

theme_update(text = element_text(size=20))
png(filename = 'coveragevsgenedensitychr1-3.png', width = 1000, height = 1000)
ggarrange(plot1, plot2, plot3)
dev.off()