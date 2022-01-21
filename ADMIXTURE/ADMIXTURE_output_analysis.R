

#copy - pasted values for CV error depending on number of clusters
K <- c(0.37549,0.36878,0.36176,0.36473,0.36300,0.36300,0.36159,0.35896,0.35685,0.35564,0.35471)

png(filename = "ADMIXTURE_cluster_number.png")
plot(K,xlab="Number of Clusters", ylab="Cross-Validation error", main="ADMIXTURE CV error by number of clusters")
dev.off()


tbl=read.table("../scripts/ADMIXTURE/diploid.3.Q")

png(filename = "ADMIXTURE_3_clusters.png", width = 2000, height = 800, units = "px")
barplot(t(as.matrix(tbl)), col=rainbow(3),
          xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()

tbl=read.table("../scripts/ADMIXTURE/diploid.10.Q")

png(filename = "ADMIXTURE_10_clusters.png", width = 2000, height = 800, units = "px")
barplot(t(as.matrix(tbl)), col=rainbow(10),
        xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()