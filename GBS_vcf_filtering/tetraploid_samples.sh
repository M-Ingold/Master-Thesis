#!/bin/bash

IN_VCF_FILE=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs_1_read_het_biallelic_SNPs_blanked_depth.vcf
OUT_VCF_FILE=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs_1_read_het_biallelic_SNPs_blanked_depth_tetra.vcf

#gets all the sample names from the vcf file, diploids deleted by hand
#grep "#CHROM" -m 1 $IN_VCF_FILE | cut -f 10- | tr "\t" "\n" > samples.txt

/home/markus/genetools/bcftools-1.14/bcftools view --samples-file  tetraploid_samples.txt --force-samples $IN_VCF_FILE > $OUT_VCF_FILE