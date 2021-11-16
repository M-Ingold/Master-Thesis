#Average genes and coverage per window tested for correlation
library(ape)

#CDS interpreted as genes
gff <- read.gff("References/DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3")
cds <- subset(gff, type == "CDS")
cds <- cds[c(1,4,5)]
cds_sorted <- cds[with(cds, order(seqid, start)), ]
cds_sorted$n <- 1
cds_chr01 <- subset(cds_sorted, seqid == "chr01")

coverage_chr01 <- read.delim("data/alignment/coverage_chr01.txt", header = F)
coverage_chr02 <- read.delim("data/alignment/coverage_chr02.txt", header = F)
coverage_chr03 <- read.delim("data/alignment/coverage_chr03.txt", header = F)

#display gene density
#hist(subset(cds_sorted, seqid == "chr01")$start)

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
#only chr1-9 so far
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
  return(df)
}

chr01genes <- rollingsum1chr(cds_chr01, 1000000)
allgenes <- rollingsumallchr(cds_sorted, 1000000)

#average coverage of the genome 
avgwindow1chr <- function(file, window){
  n = 0
  i = 0
  df <- data.frame(start=integer(), end=integer(), avg=integer())
  while (i < max(file$V2)) {
    avg <- sum(file[file$V2 %in% seq(from=i, to=i+window), ]$V3)/window
    df[n, ] <- c(i, i+window, avg)
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
plot(subset(allgenes, chr == "chr02")$sum, chr02coverage$avg)

abline(lm(chr02coverage$avg ~ subset(allgenes, chr == "chr02")$sum))

cor.test(subset(allgenes, chr == "chr02")$sum, chr02coverage$avg)



chr03coverage <- avgwindow1chr(coverage_chr03, 1000000)
plot(subset(allgenes, chr == "chr03")$sum, chr03coverage$avg)

abline(lm(chr03coverage$avg ~ subset(allgenes, chr == "chr03")$sum))

cor.test(subset(allgenes, chr == "chr03")$sum, chr03coverage$avg)

