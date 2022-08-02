#!/bin/bash

# subset VCF for SNPs with only one alternative or reference allele in all samples

IN_FILE=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs_1_read_het_biallelic_SNPs_blanked_depth.vcf
OUT_FILE=../../data/VCF/unique.vcf

# subset VCF using GATK
/home/markus/genetools/gatk-4.2.4.1/gatk SelectVariants -V $IN_FILE -O $OUT_FILE -select "AC<=1||AC>=1043" # 1044 alleles in all samples, 1043 finds samples unique for reference allele

# output file is generated from VCF
~/genetools/bcftools-1.14/bcftools stats --samples '-' $OUT_FILE | grep -E "^PSC|^# PSC" > unique_SNPs_per_sample.txt