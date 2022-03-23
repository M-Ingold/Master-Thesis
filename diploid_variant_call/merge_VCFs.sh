#!/bin/bash

VCF_FOLDER=../../data/diploid_VCF

touch filelist.txt
for entry in "$VCF_FOLDER"/freebayes_269_samples_*_diploid.vcf
do
  echo "$entry" >> filelist.txt
done

#first try bcftools for merge
bcftools concat $VCF_FOLDER/freebayes_269_samples_*_diploid.vcf --output $VCF_FOLDER/freebayes_269_samples_diploid.vcf