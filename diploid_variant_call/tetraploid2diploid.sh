#!/bin/bash

IN_FILE=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf
OUT_FILE=diploidized.vcf

cp $IN_FILE .

sed -i 's+1/1/1/1+1/1+g' *.vcf 
sed -i 's+0/1/1/1+0/1+g' *.vcf 
sed -i 's+0/0/1/1+0/1+g' *.vcf 
sed -i 's+0/0/0/1+0/1+g' *.vcf 
sed -i 's+0/0/0/0+0/0+g' *.vcf 

mv *.vcf $OUT_FILE