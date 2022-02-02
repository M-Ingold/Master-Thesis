#!/bin/bash

#####################################################################
#                                                                   #
# Split SNPs, MNPs, Indels and Complex variants into separate files #
#                                                                   #
#####################################################################

# In- and output for vcffilter
IN_VCF_FILE=$1
OUT_VCF_SNPs=$2
OUT_VCF_MNPs=$3
OUT_VCF_INDELs=$4
OUT_VCF_COMPLEX=$5

vcffilter -s -f "TYPE = snp" \
    $IN_VCF_FILE > $OUT_VCF_SNPs

vcffilter -s -f "TYPE = mnp" \
    $IN_VCF_FILE > $OUT_VCF_MNPs

vcffilter -s -f "TYPE = ins | TYPE = del" \
    $IN_VCF_FILE > $OUT_VCF_INDELs

vcffilter -s -f "TYPE = complex" \
    $IN_VCF_FILE > $OUT_VCF_COMPLEX