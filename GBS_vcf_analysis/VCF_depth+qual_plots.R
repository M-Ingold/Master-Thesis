# depth per SNP and per sample and quality per SNP are visualized before and after filtering

library(ggplot2)
library(ggpubr)

# mean depth per SNP
mean_depth_unfiltered <- read.delim("../data/VCF/Stats/unfiltered_mean_depth.ldepth.mean")
mean_depth_filtered <- read.delim("../data/VCF/Stats/filtered_mean_depth.ldepth.mean")

mean_unfiltered <- mean(mean_depth_unfiltered$MEAN_DEPTH, na.rm=T) 
mean_filtered <- mean(mean_depth_filtered$MEAN_DEPTH, na.rm=T) 

median_unfiltered <- median(mean_depth_unfiltered$MEAN_DEPTH, na.rm=T) 
median_filtered <- median(mean_depth_filtered$MEAN_DEPTH, na.rm=T) 

mean_depth_density_unfiltered <- ggplot(data = mean_depth_unfiltered, aes(x = MEAN_DEPTH)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  xlab("Mean depth per SNP") + 
  ylab("Density") +
  annotate("text", label = paste0("bar(x) == ",round(mean_unfiltered, digits = 1)), parse = TRUE, x = 150, y = 0.2, size = 5, hjust=0) + 
  annotate("text", label = paste0("tilde(x) == ",round(median_unfiltered, digits = 1)), parse = TRUE, x = 150, y = 0.18, size = 5, hjust=0) + 
  theme_bw() +
  # ggtitle("Unfiltered") +
  # theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0,200)

mean_depth_density_filtered <- ggplot(data = mean_depth_filtered, aes(x = MEAN_DEPTH)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  xlab("Mean depth per SNP") + 
  ylab("Density") +
  annotate("text", label = paste0("bar(x) == ",round(mean_filtered, digits = 1)), parse = TRUE, x = 600, y = 0.0018, size = 5, hjust=0) + 
  annotate("text", label = paste0("tilde(x) == ",round(median_filtered, digits = 1)), parse = TRUE, x = 600, y = 0.00162, size = 5, hjust=0) + 
  theme_bw()
  # ggtitle("Filtered") +
  # theme(plot.title = element_text(hjust = 0.5))
  #+xlim(0,200)


#ggarrange(mean_depth_density_unfiltered,mean_depth_density_filtered)

# mean depth per sample

mean_sample_unfiltered <- read.delim("../data/VCF/Stats/unfiltered_depth.idepth")
mean_sample_filtered <- read.delim("../data/VCF/Stats/filtered_depth.idepth")

mean_unfiltered2 <- mean(mean_sample_unfiltered$MEAN_DEPTH, na.rm=T) 
mean_filtered2 <- mean(mean_sample_filtered$MEAN_DEPTH, na.rm=T) 

median_unfiltered2 <- median(mean_sample_unfiltered$MEAN_DEPTH, na.rm=T) 
median_filtered2 <- median(mean_sample_filtered$MEAN_DEPTH, na.rm=T) 

mean_sample_density_unfiltered <- ggplot(data = mean_sample_unfiltered, aes(x = MEAN_DEPTH)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  xlab("Mean depth per sample") + 
  ylab("Density") +
  annotate("text", label = paste0("bar(x) == ",round(mean_unfiltered2, digits = 1)), parse = TRUE, x = 80, y = 0.0199, size = 5, hjust=0) + 
  annotate("text", label = paste0("tilde(x) == ",round(median_unfiltered2, digits = 1)), parse = TRUE, x = 80, y = 0.018, size = 5, hjust=0) + 
  # ggtitle("Unfiltered") +
  # theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() 
  #+xlim(0,200)

mean_sample_density_filtered <- ggplot(data = mean_sample_filtered, aes(x = MEAN_DEPTH)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  xlab("Mean depth per sample") + 
  ylab("Density") +
  annotate("text", label = paste0("bar(x) == ",round(mean_filtered2, digits = 1)), parse = TRUE, x = 610, y = 0.0018, size = 5, hjust=0) + 
  annotate("text", label = paste0("tilde(x) == ",round(median_filtered2, digits = 1)), parse = TRUE, x = 610, y = 0.00163, size = 5, hjust=0) + 
  # ggtitle("Filtered") +
  # theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() 

#ggarrange(mean_sample_density_unfiltered,mean_sample_density_filtered)

# mean quality per SNP

qual_unfiltered <- read.delim("../data/VCF/Stats/unfiltered_qual.lqual")
qual_filtered <- read.delim("../data/VCF/Stats/filtered_qual.lqual")

mean_unfiltered3 <- mean(qual_unfiltered$QUAL, na.rm=T) 
mean_filtered3 <- mean(qual_filtered$QUAL, na.rm=T) 

median_unfiltered3 <- median(qual_unfiltered$QUAL, na.rm=T) 
median_filtered3 <- median(qual_filtered$QUAL, na.rm=T) 

mean_qual_density_unfiltered <- ggplot(data = qual_unfiltered, aes(x = QUAL)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  xlab("Quality per SNP") + 
  ylab("Density") +
  annotate("text", label = paste0("bar(x) == ",round(mean_unfiltered3, digits = 1)), parse = TRUE, x = 120, y = 0.095, size = 5, hjust=0) + 
  annotate("text", label = paste0("tilde(x) == ",round(median_unfiltered3, digits = 1)), parse = TRUE, x = 120, y = 0.082, size = 5, hjust=0) + 
  theme_bw() +
  ggtitle("Unfiltered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0,200)

mean_qual_density_filtered <- ggplot(data = qual_filtered, aes(x = QUAL)) + 
  geom_density(fill = "grey", color = "black", alpha = 0.3) + 
  xlab("Quality per SNP") + 
  ylab("Density") +
  annotate("text", label = paste0("bar(x) == ",round(mean_filtered3, digits = 1)), parse = TRUE, x = 300000, y = 0.000018, size = 5, hjust=0) + 
  annotate("text", label = paste0("tilde(x) == ",round(median_filtered3, digits = 1)), parse = TRUE, x = 300000, y = 0.0000159, size = 5, hjust=0) + 
  theme_bw() +
  ggtitle("Filtered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0,500000)

#ggarrange(mean_qual_density_unfiltered,mean_qual_density_filtered)

png(filename="SNP_Depth+Qual.png", width = 750, height = 1125, res = 100)
ggarrange(mean_qual_density_unfiltered,mean_qual_density_filtered, mean_depth_density_unfiltered, mean_depth_density_filtered,  mean_sample_density_unfiltered, mean_sample_density_filtered, ncol = 2, nrow = 3)
dev.off()
