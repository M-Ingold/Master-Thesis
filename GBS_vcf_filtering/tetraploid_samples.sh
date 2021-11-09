#!/bin/bash

IN_VCF_FILE=../../data/VCF/freebayes_269_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf
OUT_VCF_FILE=../../data/VCF/freebayes_269_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_tetra.vcf

#gets all the sample names from the vcf file, diploids deleted bz hand
#grep "#CHROM" -m 1 $IN_VCF_FILE | cut -f 10- | tr "\t" "\n" > samples.txt

bcftools view --samples-file  tetraploid_samples.txt $IN_VCF_FILE > $OUT_VCF_FILE