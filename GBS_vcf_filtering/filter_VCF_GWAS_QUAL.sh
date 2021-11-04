#!/bin/bash


# In- and output for vcffilter
IN_VCF_FILE=$1
OUT_VCF_FILE=$2

# Filtering all sites with QUAL < 30
vcffilter -f " ! ( QUAL < 30 ) " \
    $IN_VCF_FILE > $OUT_VCF_FILE