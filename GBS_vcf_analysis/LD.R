library(VariantAnnotation)
library(updog)
library(ldsep)
library(sommer)
library(stringr)
library(ggplot2)
library(quantreg)
library(dplyr)
library(ggpubr)

chr2.gr <- GRanges(seqnames = "chr02", ranges=IRanges(start=1, end=100000000))
params <- ScanVcfParam(which=chr2.gr)
vcf <- readVcf("../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf.gz", param = params)

geno(header(vcf))

sizemat <- geno(vcf)$DP
refmat <- geno(vcf)$RO
ploidy <- 4

mout <- multidog(refmat = refmat, 
                 sizemat = sizemat, 
                 ploidy = ploidy, 
                 model = "norm",
                 nc = 11)

plot(mout, indices = sample(1:nrow(vcf), 3))
msub <- filter_snp(x = mout, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95)
nrow(msub$snpdf)

varnames <- paste0("logL_", 0:ploidy)
varnames

larray <- format_multidog(x = msub, varname = varnames)
dim(larray)

pmmat <- format_multidog(x = msub, varname = "postmean")

like_ld <- mldest(geno = larray, K = ploidy, type = "comp", nc=11)
plot(like_ld)

mom_ld <- mldest(geno = pmmat, K = ploidy, type = "comp", nc=11)

par(mar = c(2.4, 2.8, 0, 0) + 0.5, mgp = c(1.8, 0.6, 0))
plot(mom_ld$r2, like_ld$r2, 
     xlab = expression(paste(textstyle(Naive), ~~hat(r)^2)), 
     ylab = expression(paste(textstyle(MLE), ~~hat(r)^2)), 
     pch  = 20)
abline(0, 1, lty = 2, col = 2)

ldmat <- format_lddf(obj = like_ld, element = "LD")
ldmat[1:4, 1:4]

map <- data.frame(Locus=msub$snpdf$snp, LG="chr02", Position=as.numeric(str_match(msub$snpdf$snp, ":\\s*(.*?)\\s*_")[,2]))

#LD.decay(ldmat, map = map)

ld_df <- data.frame(i=map$Position[like_ld$i], j=map$Position[like_ld$j], r2=like_ld$r2, Dprime=like_ld$Dprime)

ld_df$delta <- abs(ld_df$i - ld_df$j)

plot(ld_df$delta, ld_df$r2)

plot(ld_df$delta, ld_df$Dprime)
# Dat.nls <- nls(r2 ~ SSlogis(delta, Asym, mid, scal), data=ld_df, control=nls.control(minFactor=1/1000000000,maxiter=100,warnOnly=T))
# ; Dat.nls
# lines(1:25, predict(Dat.nls, newdata=list(x=1:25)), col=1)

# average r^2 plotted
n<-264
distance<-ld_df$delta
LD.data<-ld_df$r2

file <- data.frame(dist = distance, rsq = LD.data)
file <- na.omit(file)
file$distc <- cut(file$dist,breaks=seq(from=min(file$dist)-1,to=max(file$dist)+1,by=100000))
dfr1 <- file %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
gg1 <- ggplot()+
  geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8")) +
  ggtitle("r^2 mean")

gg2 <- ggplot()+
  geom_point(data=dfr1,aes(x=start,y=median),size=0.4,colour="grey20")+
  geom_line(data=dfr1,aes(x=start,y=median),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8")) +
  ggtitle("r^2 median")

theme_update(text = element_text(size=20))
png(filename = 'rsq_raw.png', width = 1000, height = 500)
ggarrange(gg1, gg2)
dev.off()

#file <- subset(file, rsq > 0.01)
#formula <- ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance)))

# HW.st<-c(C=0.1)
# HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=1000))
# tt<-summary(HW.nonlinear)
# new.rho<-tt$parameters[1]
# fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
# 
# model_result <- summary(HW.nonlinear)
# cor(ld_df$r2, predict(HW.nonlinear))
# 
# lines(ld_df$delta, order(fpoints))
# 
# HW.nonlinear <- nlrq()

# C values range from about 0.5 to 2, start with 0.1
Cstart <- c(C=0.1)

# fit a non linear model using the arbitrary C value, 
# N is the number of the genotypes that have the SNP site
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*n*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*n*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile <- data.frame(file$dist, newrsq)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]

# plotting the values
png("LD_decay.png", height=500, width = 500)
#mar.default <- c(5,4,4,2) + 0.1
par(mar = c(5, 5, 1, 1)) 
plot(file$dist, file$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile$file.dist, newfile$newrsq, col="red", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizental line
abline(v=halfdecaydist, col="green")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=1, col="green")
dev.off()

# nlrq 
fit.nlrq <- nlrq(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                    ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*n*(2+C*dist) * (11+C*dist) ) ), 
                  data=file, start = Cstart, tau = 0.9)

summary(fit.nlrq)

# Doing the same thing as above leads to single-digit LD values
# rho <- -0.9
# 
# newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
#   ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
#       (2*n*(2+rho*file$dist) * (11+rho*file$dist) ) )
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

plot(file$dist, file$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")

lines(file$dist, predict(fit.nlrq), lty = 2, col = "red")
      