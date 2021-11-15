#Rolling gene average from GFF file
library(ape)
library(zoo)

gff <- read.gff("References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3")

cds <- subset(gff, type == "CDS")
cds <- cds[c(1,4,5)]

cds_sorted <- cds[with(cds, order(seqid, start)), ]
cds_sorted$n <- 1

cds_chr01 <- subset(cds_sorted, seqid == "chr01")
#hist(subset(cds_sorted, seqid == "chr01")$start)

#rolling average for coverage, here count per window size?

rollingsum <- function(file, window){
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

rollingsum(cds_chr01, 1000000)
# 
# n = 0
# i = 0
# window = 1000000
# while (i < max(cds_chr01$end)) {
#   sum <- sum(cds_chr01[cds_chr01$start %in% seq(from=i, to=i+window), ]$n)
#     #sum(cds_chr01[cds_chr01$pos > i || cds_chr01$pos < i+window, ]$n)
#   df[n, ] <- c(i, i+window, sum)
#   n <- n+1
#   i <- i+window
# }
# 
# sum(cds_chr01[cds_chr01$start %in% seq(from=i, to=i+window), ]$n)
