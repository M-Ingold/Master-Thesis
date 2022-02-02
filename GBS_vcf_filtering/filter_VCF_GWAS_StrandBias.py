#!/bin/python

# importing modules
import pathlib
import vcf
import pandas as pd
import numpy as np
import sys

########################
#                      #
# Reading the VCF file #
#                      #
########################

# Relative input string for the pyVCF module
relVcfPathString = sys.argv[1]
# Turning it into absolute path
vcfPathString = pathlib.Path(relVcfPathString)
absVcfPath = vcfPathString.resolve(strict=True)
absVcfPath = str(absVcfPath)
# Reading it as VCF file
vcfFile = vcf.Reader(filename=absVcfPath)

############################
#                          #
# Creating output VCF file #
#                          #
############################

# Relative input string for the pyVCF module
relVcfOutPathString = sys.argv[2]
# Turning it into absolute path
vcfOutPathString = pathlib.Path(relVcfOutPathString)
absVcfOutPath = vcfOutPathString.resolve(strict=True)
absVcfOutPath = str(absVcfOutPath)
# Creating the output object using the output string
vcfOutFile = vcf.Writer(open(absVcfOutPath, 'w'), vcfFile)

#############################################################################################
#                                                                                           #
# Actual filtering and writing variants that passed the missing sample percentage threshold #
#                                                                                           #
#############################################################################################
"""
Variants were the alleles were only observed on either the forward or the reverse strand 
can be subjected to strand bias which is a kind of sequencing error, where 
"""

# Filtering the SNPs according to the minimum number of reference and alternative allele
for rec in vcfFile:
    # Check if the reference allele was observed at this position
    if rec.INFO['RO'] > 0:
        # If both, REF and ALT allele are observed at least once on both forward and reverse strand, no strand bias is assumed and the SNP passes the filter
        if (rec.INFO['SRF'] > 0 and rec.INFO['SRR'] > 0) and (rec.INFO['SAF'][0] > 0 and rec.INFO['SAR'][0] > 0):
            vcfOutFile.write_record(rec)
    # Check if there are no REF observations but two ALT alleles observed at this position
    #elif rec.INFO['RO'] == 0 and len(rec.INFO["AO"]) == 2:
        # If both ALT alleles are observed at least once on both forward and reverse strand, no strand bias is assumed and the SNP passes the filter
        #if (rec.INFO['SAF'][0] > 0 and rec.INFO['SAR'][0] > 0) and (rec.INFO['SAF'][1] > 0 and rec.INFO['SAR'][1] > 0):
        #    vcfOutFile.write_record(rec)