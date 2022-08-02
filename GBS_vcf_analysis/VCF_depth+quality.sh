#!/bin/bash

# extract overall depth, mean depth per sample and quality from unfiltered and filtered vcf for plotting in R

UNFILTERED_VCF=../../data/VCF/freebayes_269_samples_chr01-12.vcf
FILTERED_VCF=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs_1_read_het_biallelic_SNPs_blanked_depth.vcf

#create folter for stats
TARGET_DIR=../../data/VCF/Stats
mkdir -p $TARGET_DIR

# extract desired data with vcftools and save them in specified folder
vcftools --vcf $UNFILTERED_VCF --site-mean-depth --out $TARGET_DIR/unfiltered_mean_depth
vcftools --vcf $UNFILTERED_VCF --depth --out $TARGET_DIR/unfiltered_depth
vcftools --vcf $UNFILTERED_VCF --site-quality --out $TARGET_DIR/unfiltered_qual

vcftools --vcf $FILTERED_VCF --site-mean-depth --out $TARGET_DIR/filtered_mean_depth
vcftools --vcf $FILTERED_VCF --depth --out $TARGET_DIR/filtered_depth
vcftools --vcf $FILTERED_VCF --site-quality --out $TARGET_DIR/filtered_qual
