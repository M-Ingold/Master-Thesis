#!/bin/bash

IN_FILE=../../data/diploid_VCF/diploidized.vcf
OUT_FILE=../../data/diploid_VCF/structure

~/genetools/plink_linux_x86_64_20210606/plink --vcf $IN_FILE --recode structure --out $OUT_FILE