ld <- read.delim("Lddecay.stat",stringsAsFactors = FALSE,header=TRUE, sep = "\t") 

file <- ld

colnames(file) <- c("dist","rsq","D","Sumrsq","SumD","numofpairs")

# C values range from about 0.5 to 2, start with 0.1
Cstart <- c(C=0.1)
N<-261
# fit a non linear model using the arbitrary C value, 
# N is the number of the genotypes that have the SNP site
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))

# modelC2 <- nlrq(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
#                   ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
#                 data=file, start=Cstart, tau = 0.1, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
### extremely low! ###
rho <- summary(modelC)$parameters[1]
#rho2 <- summary(modelC2)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile <- data.frame(file$dist, newrsq)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]

# plotting the values
pdf("LD_decay.pdf", height=5, width = 5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file$dist, file$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile$file.dist, newfile$newrsq, col="red", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizental line
abline(v=halfdecaydist, col="green")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
dev.off()
