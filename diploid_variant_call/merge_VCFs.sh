#!/bin/bash

VCF_FOLDER=../../data/diploid_VCF

touch filelist.txt
for entry in "$VCF_FOLDER"/*
do
  echo "$VCF_FOLDER/$entry" > filelist.txt
done

#first try bcftools for merge
bcftools concat filelist.txt --output freebayes_269_samples_${chr}_diploid.vcf