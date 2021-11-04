#!/bin/bash

# In- and output file paths
IN_PATH=$1
OUT_PATH=$2

# Filtering for allele bias of alternative alleles
vcffilter -s -f "SAF > 0 & SAR > 0" $IN_PATH > $OUT_PATH