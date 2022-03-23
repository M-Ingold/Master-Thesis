library(VariantAnnotation)
library(updog)
library(ldsep)
library(sommer)
library(stringr)
library(ggplot2)
library(quantreg)
library(dplyr)
library(ggpubr)

# chr2.gr <- GRanges(seqnames = "chr02", ranges=IRanges(start=1, end=100000000))
# params <- ScanVcfParam(which=chr2.gr)
# vcf <- readVcf("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf.gz", param = params)

# create list of chromosomes
ploidy <- 4

chr <- c()
for (i in 1:9) {
  chr <- c(chr, paste0("chr0", i))
}
for (i in 10:12) {
  chr <- c(chr, paste0("chr", i))
}

# create GRanges + parameters for reading in data by chromosome
param_list <- c()
for (i in chr) {
  param_list <- c(param_list, ScanVcfParam(which=GRanges(seqnames = i, ranges=IRanges(start=1, end=100000000))))
}

# create list of SNPs by chromosome
vcf_list_by_chr <- list()
n <- 1
for (param in param_list) {
  vcf_list_by_chr[[n]] <- readVcf("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf.gz", param = param)
  n <- n+1
}


# run updog for each chromosome and filter out mostly monoallelic SNPs (synonymus to filtering MAF?)
multidog_results <- list()
for (i in 1:12) {
  sizemat <- geno(vcf_list_by_chr[[i]])$DP
  refmat <- geno(vcf_list_by_chr[[i]])$RO

  x <- multidog(refmat = refmat, 
                sizemat = sizemat, 
                ploidy = 4, 
                model = "norm",
                nc = 11)
  msub <- filter_snp(x = x, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95)
  multidog_results[[i]] <- msub
}

# format multidog results and input them into LD estimator function
varnames <- paste0("logL_", 0:ploidy)
like_lds <- list()
for (i in 1:12) {
  larray <- format_multidog(x = multidog_results[[i]], varname = varnames)
  like_lds[[i]] <- mldest(geno = larray, K = ploidy, type = "comp", nc=11, se = F) # se = FALSE for shorter computation time
}

# repeat for chr 12, delete later
larray <- format_multidog(x = multidog_results[[12]], varname = varnames)
like_lds[[12]] <- mldest(geno = larray, K = ploidy, type = "comp", nc=11, se = F)




# save like_lds in case of crash
save(like_lds, file = "like_lds")

# create r2 matrix
ldmats <- list()
for (i in 1:12) {
  ldmats[[i]] <- format_lddf(obj = like_lds[[i]], element = "r2")
}

# extract SNP positions from multidog_results
maps <- list()
for (i in 1:12) {
  maps[[i]] <- data.frame(Locus=multidog_results[[i]]$snpdf$snp, Position=as.numeric(str_match(multidog_results[[i]]$snpdf$snp, ":\\s*(.*?)\\s*_")[,2]))
}

# create LD dataframe
ld_dfs <- list()
for (n in 1:12) {
  df <- data.frame(i=maps[[n]]$Position[like_lds[[n]]$i], j=maps[[n]]$Position[like_lds[[n]]$j], r2=like_lds[[n]]$r2, Dprime=like_lds[[n]]$Dprime)

  df$dist <- abs(df$i - df$j)
  df <- df[order(df$dist),]
  ld_dfs[[n]] <- df
}

# Method according to Hill WG, Weir BS (1988)
# estimate 90th percentile C of the LD decay formula

# Start value for estimating C
Cstart <- c(C=0.01)

newfiles <- list()
for (i in 1:12) {
  fit.nlrq <- nlrq(r2 ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( n*(2+C*dist) * (11+C*dist) ) ),
                 data=ld_dfs[[i]], start = Cstart, tau = .9, nlrq.control(InitialStepSize = 0.01, maxiter = 100), trace = F)
  
  # extract the recombination parameter in 4Nr
  rec <- summary(fit.nlrq)$coefficients[1]
  
  newrsq <- ( (10+rec*ld_dfs[[i]]$dist) / ( (2+rec*ld_dfs[[i]]$dist) * (11+rec*ld_dfs[[i]]$dist) ) ) *
    ( 1 + ( (3+rec * ld_dfs[[i]]$dist) * (12+12*rec*ld_dfs[[i]]$dist + (rec*ld_dfs[[i]]$dist)^2) ) / 
    (n*(2+rec*ld_dfs[[i]]$dist) * (11+rec*ld_dfs[[i]]$dist) ) )
  nf <- data.frame(ld_dfs[[i]]$dist, newrsq)
  newfiles[[i]] <- nf[order(nf$ld_dfs..i...dist),]
}

# plot LD
par(mfrow=c(3,4))
halfdecaydist <- list()
tenthdecaydist <- list()
for (i in 1:12) {
  maxld <- max(newfiles[[i]]$newrsq,na.rm=TRUE) #using max LD value from adjusted data
  halfdecay = 0.5 * maxld
  halfdecaydist[[i]] <- newfiles[[i]]$ld_dfs..i...dist[which.min(abs(newfiles[[i]]$newrsq-halfdecay))]
  tenthdecay = 0.1
  tenthdecaydist[[i]] <- newfiles[[i]]$ld_dfs..i...dist[which.min(abs(newfiles[[i]]$newrsq-tenthdecay))]
  
  # plotting the values
  par(mar = c(5, 5, 1, 1)) 
  plot(ld_dfs[[i]]$dist, ld_dfs[[i]]$r2, pch=".", cex=2, xlab="Distance (bp)", main=paste0("Chromosome",i), ylab=expression(LD ~ (r^2)), col="grey")
  lines(newfiles[[i]]$ld_dfs..i...dist, newfiles[[i]]$newrsq, col="red", lwd=2)
  abline(h=0.1, col="blue") # if you need to add horizental line
  abline(v=tenthdecaydist[[i]], col="green")
  mtext(round(tenthdecaydist[[i]],2), side=1, line=0.05, at=tenthdecaydist[[i]], cex=1, col="green")
  #lines(ld_df$dist, predict(spline_fit,), lty = 2, col = "black") # add spline
}


# sizemat <- geno(vcf)$DP
# refmat <- geno(vcf)$RO
# ploidy <- 4
# 
# mout <- multidog(refmat = refmat,
#                  sizemat = sizemat,
#                  ploidy = ploidy,
#                  model = "norm",
#                  nc = 11)
# 
# plot(mout, indices = sample(1:nrow(vcf), 3))
# msub <- filter_snp(x = mout, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95)
# nrow(msub$snpdf)

# varnames <- paste0("logL_", 0:ploidy)
# varnames
# 
# larray <- format_multidog(x = msub, varname = varnames)
# dim(larray)
# 
# pmmat <- format_multidog(x = msub, varname = "postmean")
# 
# like_ld <- mldest(geno = larray, K = ploidy, type = "comp", nc=11)
# plot(like_ld)
# 
# mom_ld <- mldest(geno = pmmat, K = ploidy, type = "comp", nc=11)
# 
# par(mar = c(2.4, 2.8, 0, 0) + 0.5, mgp = c(1.8, 0.6, 0))
# plot(mom_ld$r2, like_ld$r2, 
#      xlab = expression(paste(textstyle(Naive), ~~hat(r)^2)), 
#      ylab = expression(paste(textstyle(MLE), ~~hat(r)^2)), 
#      pch  = 20)
# abline(0, 1, lty = 2, col = 2)

# ldmat <- format_lddf(obj = like_ld, element = "r2")
# ldmat[1:4, 1:4]

# map <- data.frame(Locus=msub$snpdf$snp, LG="chr02", Position=as.numeric(str_match(msub$snpdf$snp, ":\\s*(.*?)\\s*_")[,2]))

##LD.decay(ldmat, map = map)

# ld_df <- data.frame(i=map$Position[like_ld$i], j=map$Position[like_ld$j], r2=like_ld$r2, Dprime=like_ld$Dprime)
# 
# ld_df$dist <- abs(ld_df$i - ld_df$j)

#plot(ld_df$dist, ld_df$r2)

#plot(ld_df$dist, ld_df$Dprime)

# SSlogis doesn't work (value becomes too small?)
#Dat.nls <- nls(r2 ~ SSlogis(delta, Asym, mid, scal), data=ld_df, control=nls.control(minFactor=1/1000000000,maxiter=100,warnOnly=T))
# ; Dat.nls
# lines(1:25, predict(Dat.nls, newdata=list(x=1:25)), col=1)

## average r^2 plotted
# n<-261
# distance<-ld_df$delta
# LD.data<-ld_df$r2
# 
# file <- data.frame(dist = distance, rsq = LD.data)
# file <- na.omit(file)
# file$distc <- cut(file$dist,breaks=seq(from=min(file$dist)-1,to=max(file$dist)+1,by=100000))
# dfr1 <- file %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
# dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
#                         end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
#                         mid=start+((end-start)/2))
# gg1 <- ggplot()+
#   geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
#   geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
#   labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
#   scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8")) +
#   ggtitle("r^2 mean")
# 
# gg2 <- ggplot()+
#   geom_point(data=dfr1,aes(x=start,y=median),size=0.4,colour="grey20")+
#   geom_line(data=dfr1,aes(x=start,y=median),size=0.3,alpha=0.5,colour="grey40")+
#   labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
#   scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8")) +
#   ggtitle("r^2 median")
# 
# theme_update(text = element_text(size=20))
# png(filename = 'rsq_raw.png', width = 1000, height = 500)
# ggarrange(gg1, gg2)
# dev.off()



# Method according to Hill WG, Weir BS (1988)
# pretty sure it was done differently by Vos and Sharma, because the y-Axis intercept of the function below can't be > 0.5

# Start value for estimating C
#Cstart <- c(C=0.01)

# fit a non linear model using the arbitrary C value, 
# n is the number of the genotypes that have the SNP site
# modelC <- nls(r2 ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
#                 ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( n*(2+C*dist) * (11+C*dist) ) ), 
#               data=ld_df, start=Cstart, control=nls.control(maxiter=100), trace = T)

# nlrq 
#ld_df <- ld_df[order(ld_df$dist),]

# fit.nlrq <- nlrq(r2 ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( n*(2+C*dist) * (11+C*dist) ) ),
#                   data=ld_df, start = Cstart, tau = .9, nlrq.control(InitialStepSize = 0.01, maxiter = 100), trace = T)
# 
# summary <- summary(fit.nlrq)
# #summary$parameters[1]
# summary$coefficients[1]
# # extract the recombination parameter in 4Nr
# rec <- summary(modelC)$parameters[1]
# rec <- 3.31994e-06 # why did I get this earlier?
# rec <- 1.188454e-05 
# 
# # feed in the new value of rec to obtain LD values adjusted for their distances along the chromosome/genome
# # shouldn't this be the fit of the nlrq function?
# newrsq <- ( (10+rec*ld_df$dist) / ( (2+rec*ld_df$dist) * (11+rec*ld_df$dist) ) ) *
#   ( 1 + ( (3+rec * ld_df$dist) * (12+12*rec*ld_df$dist + (rec*ld_df$dist)^2) ) / 
#       (n*(2+rec*ld_df$dist) * (11+rec*ld_df$dist) ) )

#newfile <- data.frame(ld_df$dist, newrsq)

# maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
# halfdecay = 0.5 * maxld
# 
# halfdecaydist <- newfile$ld_df.dist[which.min(abs(newfile$newrsq-halfdecay))]
# newfile <- newfile[order(newfile$ld_df.dist),]
# # add 1/10 decay
# tenthdecay = 0.1
# tenthdecaydist <- newfile$ld_df.dist[which.min(abs(newfile$newrsq-tenthdecay))]
# 
# # plotting the values
# png("LD_decay.png", height=500, width = 500)
# #mar.default <- c(5,4,4,2) + 0.1
# par(mar = c(5, 5, 1, 1)) 
# plot(ld_df$dist, ld_df$r2, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
# lines(newfile$ld_df.dist, newfile$newrsq, col="red", lwd=2)
# abline(h=0.1, col="blue") # if you need to add horizental line
# abline(v=tenthdecaydist, col="green")
# mtext(round(tenthdecaydist,2), side=1, line=0.05, at=tenthdecaydist, cex=1, col="green")
# #lines(ld_df$dist, predict(spline_fit,), lty = 2, col = "black") # add spline
# dev.off()



# Doing the same thing as above leads to single-digit LD values
# rec <- -0.9
# 
# newrsq <- ( (10+rec*file$dist) / ( (2+rec*file$dist) * (11+rec*file$dist) ) ) *
#   ( 1 + ( (3+rec * file$dist) * (12+12*rec*file$dist + (rec*file$dist)^2) ) / 
#       (2*n*(2+rec*file$dist) * (11+rec*file$dist) ) )
# 
# newfile <- data.frame(file$dist, newrsq)
# 
# maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
# halfdecay = maxld*0.5
# halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
# newfile <- newfile[order(newfile$file.dist),]
# 
# 
# plot(file$dist, file$rsq, cex=2, pch=".", xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
# 
# 
# par(mar = c(5, 5, 1, 1)) 
# plot(file$dist, file$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
# lines(newfile$file.dist, newfile$newrsq, col="red", lwd=2)
# abline(h=0.1, col="blue") # if you need to add horizental line
# abline(v=halfdecaydist, col="green")
# mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=1, col="green")

#x <- seq(min(file$dist), max(file$dist), length = 10000)

# plot(ld_df$dist, ld_df$r2, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
# lines(ld_df$dist, predict(modelC,), lty = 2, col = "blue")
# 
# lines(ld_df$dist, predict(fit.nlrq,), lty = 2, col = "red")
      

library(splines)
# spline 90th percentile regression, pretty sure this is not the right way to do this
# try different degree of freedoms?
spline_fit <- rq(r2 ~ bs(dist, df = 5), tau = 0.9, data = ld_df)
lines(ld_df$dist, predict(spline_fit,), lty = 2, col = "black")
# 
# library(gaston)
# # LD heatmap, doesn't work for this large dataset
# #library(Matrix)
# #like_ld_sym <- as.matrix(forceSymmetric(ldmat))
# png(filename = "test.png", width = 1000, height = 1000)
# LD.plot(ldmat, map$Position, write.snp.id = F, write.ld = NULL, draw.chr = F
#         , pdf.file = "test.pdf", finalize.pdf = T
#         )
# dev.off()

#library(LDheatmap) # package ‘LDheatmap’ is not available (for R version 3.6.3)

### no heatmap for now ###
