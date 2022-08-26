# LD and LD decay was determined and plotted using this script

library(VariantAnnotation)
library(updog)
library(ldsep)
library(sommer)
library(stringr)
library(ggplot2)
library(quantreg)
library(dplyr)
library(ggpubr)
library(splines)

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

# plot three random SNP's multidog results. A chromosome number needs to be specified in the brackets
plot(multidog_results[[2]], indices = sample(1:nrow(vcf_list_by_chr[[2]]), 3))
# a nice example plot was saved using export

# format multidog results and input them into LD estimator function
varnames <- paste0("logL_", 0:ploidy)
like_lds <- list()
for (i in 1:12) {
  larray <- format_multidog(x = multidog_results[[i]], varname = varnames)
  like_lds[[i]] <- mldest(geno = larray, K = ploidy, type = "comp", nc=11, se = F) # se = FALSE for shorter computation time
}



# save like_lds in case of crash
save(like_lds, file = "like_lds")
#load("like_lds")

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
  df <- na.omit(df[order(df$dist),])
  ld_dfs[[n]] <- df
}

# Method according to Hill WG, Weir BS (1988), this wasn't used for the master thesis but was the first attempt at estimating LD decay
# estimate 90th percentile C of the LD decay formula

# Start value for estimating C
# Cstart <- c(C=0.01)
# 
# newfiles <- list()
# for (i in 1:12) {
#   fit.nlrq <- nlrq(r2 ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( n*(2+C*dist) * (11+C*dist) ) ),
#                  data=ld_dfs[[i]], start = Cstart, tau = .9, nlrq.control(InitialStepSize = 0.01, maxiter = 100), trace = F)
#   
#   # extract the recombination parameter in 4Nr
#   rec <- summary(fit.nlrq)$coefficients[1]
#   
#   newrsq <- ( (10+rec*ld_dfs[[i]]$dist) / ( (2+rec*ld_dfs[[i]]$dist) * (11+rec*ld_dfs[[i]]$dist) ) ) *
#     ( 1 + ( (3+rec * ld_dfs[[i]]$dist) * (12+12*rec*ld_dfs[[i]]$dist + (rec*ld_dfs[[i]]$dist)^2) ) / 
#     (n*(2+rec*ld_dfs[[i]]$dist) * (11+rec*ld_dfs[[i]]$dist) ) )
#   nf <- data.frame(ld_dfs[[i]]$dist, newrsq)
#   newfiles[[i]] <- nf[order(nf$ld_dfs..i...dist),]
# }

# # plot LD from formula
# par(mfrow=c(3,4))
# halfdecaydist <- list()
# tenthdecaydist <- list()
# for (i in 1:12) {
#   maxld <- max(newfiles[[i]]$newrsq,na.rm=TRUE) #using max LD value from adjusted data
#   halfdecay = 0.5 * maxld
#   halfdecaydist[[i]] <- newfiles[[i]]$ld_dfs..i...dist[which.min(abs(newfiles[[i]]$newrsq-halfdecay))]
#   tenthdecay = 0.1
#   tenthdecaydist[[i]] <- newfiles[[i]]$ld_dfs..i...dist[which.min(abs(newfiles[[i]]$newrsq-tenthdecay))]
#   
#   # plotting the values
#   par(mar = c(5, 5, 1, 1)) 
#   plot(ld_dfs[[i]]$dist, ld_dfs[[i]]$r2, pch=".", cex=2, xlab="Distance (bp)", main=paste0("Chromosome",i), ylab=expression(LD ~ (r^2)), col="grey")
#   lines(newfiles[[i]]$ld_dfs..i...dist, newfiles[[i]]$newrsq, col="red", lwd=2)
#   abline(h=0.1, col="blue") # if you need to add horizental line
#   abline(v=tenthdecaydist[[i]], col="green")
#   mtext(round(tenthdecaydist[[i]],2), side=1, line=0.05, at=tenthdecaydist[[i]], cex=1, col="green")
#   #lines(ld_df$dist, predict(spline_fit,), lty = 2, col = "black") # add spline
# }

# testing for good df values. Plot with desired spline fits was saved manually
spline_fit3 <- rq(r2 ~ bs(dist, df=3), tau = 0.9, data = ld_df)
spline_fit4 <- rq(r2 ~ bs(dist, df=4), tau = 0.9, data = ld_df)
spline_fit5 <- rq(r2 ~ bs(dist, df=5), tau = 0.9, data = ld_df)
spline_fit6 <- rq(r2 ~ bs(dist, df=6), tau = 0.9, data = ld_df)
spline_fit7 <- rq(r2 ~ bs(dist, df=7), tau = 0.9, data = ld_df)
spline_fit8 <- rq(r2 ~ bs(dist, df=8), tau = 0.9, data = ld_df)
spline_fit9 <- rq(r2 ~ bs(dist, df=9), tau = 0.9, data = ld_df)
spline_fit10 <- rq(r2 ~ bs(dist, df=10), tau = 0.9, data = ld_df)

spline_fit20 <- rq(r2 ~ bs(dist, df=20), tau = 0.9, data = ld_df)

spline_fit30 <- rq(r2 ~ bs(dist, df=30), tau = 0.9, data = ld_df)

plot(ld_df$dist, ld_df$r2, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey", main = "90 percentile splines with different DFs")
lines(ld_df$dist, predict(spline_fit3,), lty = 1, col = "black", lwd = 2)
lines(ld_df$dist, predict(spline_fit4,), lty = 1, col = "green", lwd = 2)
lines(ld_df$dist, predict(spline_fit5,), lty = 1, col = "red", lwd = 2)
lines(ld_df$dist, predict(spline_fit6,), lty = 1, col = "blue", lwd = 2)
lines(ld_df$dist, predict(spline_fit7,), lty = 2, col = "green")
lines(ld_df$dist, predict(spline_fit8,), lty = 2, col = "yellow")
lines(ld_df$dist, predict(spline_fit9,), lty = 2, col = "purple")
lines(ld_df$dist, predict(spline_fit10,), lty = 2, col = "black")

lines(ld_df$dist, predict(spline_fit20,), lty = 1, col = "orange")
lines(ld_df$dist, predict(spline_fit30,), lty = 2, col = "cyan")



# spline fit as in Vos and Sharma's publication. As no parameter was given in other papers, df = 5 was guessed based on the curve stabilizing at high distances compared to other values
spline_fits <- list()
for (i in 1:12) {
  spline_fit <- rq(r2 ~ bs(dist, df = 5), tau = 0.9, data = ld_dfs[[i]])
  spline_fits[[i]] <- data.frame(ld_dfs[[i]]$dist, predict(spline_fit))
}

# plot LD from spline
par(mfrow=c(6,4))
halfdecaydist_spline <- list()
tenthdecaydist_spline <- list()
for (i in 1:12) {
  maxld <- max(spline_fits[[i]][2],na.rm=TRUE) 
  halfdecay = 0.5 * maxld
  halfdecaydist_spline[[i]] <- spline_fits[[i]][subset(spline_fits[[i]][2])<halfdecay,][1,1]
  tenthdecay = 0.1
  tenthdecaydist_spline[[i]] <- spline_fits[[i]][subset(spline_fits[[i]][2])<tenthdecay,][1,1]
  
  # plotting the values
  par(mar = c(5, 5, 1, 1)) 
  plot(ld_dfs[[i]]$dist, ld_dfs[[i]]$r2, pch=".", cex=2, xlab="Distance (bp)", main=paste0("Chromosome",i), ylab=expression(LD ~ (r^2)), col="grey")
  lines(spline_fits[[i]], lwd = 2, col = "red")
  abline(h=0.1, col="blue") # if you need to add horizental line
  #abline(v=tenthdecaydist_spline[[i]], col="green") # tenthdecaydistance line
  #mtext(round(tenthdecaydist_spline[[i]],2), side=1, line=0.05, at=tenthdecaydist_spline[[i]], cex=1, col="green") # tenthdecaydistance number
  #lines(ld_df$dist, predict(spline_fit,), lty = 2, col = "black") # add spline
}
# plot saved manually

# create dataframe with all decay distances
decay_distances <- data.frame(row.names = chr, half = round(unlist(halfdecaydist_spline)/10^6, digits = 2), tenth = round(unlist(tenthdecaydist_spline)/10^6, digits = 2))
# output as latex table
library(xtable)
xtable(decay_distances, auto = T)


# plot the comparison of rsq based on variant call genotypes vs. based on genotype probabilities
pmmat <- format_multidog(x = msub, varname = "postmean")

mom_ld <- mldest(geno = pmmat, K = ploidy, type = "comp", nc=11)

par(mar = c(2.4, 2.8, 0, 0) + 0.5, mgp = c(1.8, 0.6, 0))
plot(mom_ld$r2, like_ld$r2,
     xlab = expression(paste(textstyle(Naive), ~~ r^2)),
     ylab = expression(paste(textstyle(Likelihood), ~~ r^2)),
     pch  = 20)
abline(0, 1, lty = 2, col = 2)
# plot saved manually



# plot mean and median rsq per distance. Not used in master thesis
# ld_df <- data.frame(i=map$Position[like_ld$i], j=map$Position[like_ld$j], r2=like_ld$r2, Dprime=like_ld$Dprime)
# 
# ld_df$dist <- abs(ld_df$i - ld_df$j)

#plot(ld_df$dist, ld_df$r2)

#plot(ld_df$dist, ld_df$Dprime)


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
# png(filename = 'rsq_raw.png', width = 1000, height = 500)
# ggarrange(gg1, gg2)
# dev.off()


 

# library(gaston)
# gaston LD heatmap, doesn't work for this large dataset
# 
# png(filename = "test.png", width = 1000, height = 1000)
# LD.plot(like_lds[[1]], maps[[1]]$Position, write.snp.id = F, write.ld = NULL, draw.chr = T
#         , pdf.file = "test.pdf", finalize.pdf = T, max.dist = F
#         )
# dev.off()

library(LDheatmap) 
library(grid)
library(gridExtra)
library(lattice)

rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")

#png(filename = "test.png", width = 2000, height = 2000)
#par(mfrow=c(3,4)) # doesn't work with this

# LDheatmaps can't be stored in lists???
# plot_list <- list()
# for (i in 1:12) {
#   plot_list[i] <-LDheatmap(ldmats[[i]],maps[[i]]$Position, flip = T, title = paste("Chromosome", i),text = F, color=rgb.palette(30), pop = F)$LDheatmapGrob
#   }
# grid.arrange(plot_list, ncol = 4)




LDheat1 <-LDheatmap(ldmats[[1]],maps[[1]]$Position, flip = F, title = paste("Chromosome", 1),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat1 <- editGrob(LDheat1$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat2 <-LDheatmap(ldmats[[2]],maps[[2]]$Position, flip = F, title = paste("Chromosome", 2),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat2 <- editGrob(LDheat2$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat3 <-LDheatmap(ldmats[[3]],maps[[3]]$Position, flip = F, title = paste("Chromosome", 3),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat3 <- editGrob(LDheat3$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat4 <-LDheatmap(ldmats[[4]],maps[[4]]$Position, flip = F, title = paste("Chromosome", 4),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat4 <- editGrob(LDheat4$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat5 <-LDheatmap(ldmats[[5]],maps[[5]]$Position, flip = F, title = paste("Chromosome", 5),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat5 <- editGrob(LDheat5$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat6 <-LDheatmap(ldmats[[6]],maps[[6]]$Position, flip = F, title = paste("Chromosome", 6),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat6 <- editGrob(LDheat6$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat7 <-LDheatmap(ldmats[[7]],maps[[7]]$Position, flip = F, title = paste("Chromosome", 7),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat7 <- editGrob(LDheat7$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat8 <-LDheatmap(ldmats[[8]],maps[[8]]$Position, flip = F, title = paste("Chromosome", 8),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat8 <- editGrob(LDheat8$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat9 <-LDheatmap(ldmats[[9]],maps[[9]]$Position, flip = F, title = paste("Chromosome", 9),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat9 <- editGrob(LDheat9$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat10 <-LDheatmap(ldmats[[10]],maps[[10]]$Position, flip = F, title =paste("Chromosome", 10),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat10 <- editGrob(LDheat10$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat11 <-LDheatmap(ldmats[[11]],maps[[11]]$Position, flip = F, title = paste("Chromosome", 11),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat11 <- editGrob(LDheat11$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))
LDheat12 <-LDheatmap(ldmats[[12]],maps[[12]]$Position, flip = F, title = paste("Chromosome", 12),text = F, color=rgb.palette(30), pop = F)#$LDheatmapGrob
LDheat12 <- editGrob(LDheat12$LDheatmapGrob, gPath("heatMap", "title"),  gp=gpar(cex=10))

# pdf(file = "physical_LD.pdf", width = 100, height = 80)
# grid.arrange(LDheat1,LDheat2,LDheat3,LDheat4,LDheat5,LDheat6,LDheat7,LDheat8,LDheat9,LDheat10,LDheat11,LDheat12, ncol =4)
# dev.off()

png(filename = "physical_LD.png", width = 6000, height = 8000)
grid.arrange(LDheat1,LDheat2,LDheat3,LDheat4,LDheat5,LDheat6,LDheat7,LDheat8,LDheat9,LDheat10,LDheat11,LDheat12, ncol =3)
dev.off()


