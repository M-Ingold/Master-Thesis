#!/bin/python

# importing modules
import pathlib
import vcf
import pandas as pd
import numpy as np
import sys

##############################################################################
#                                                                            #
# Reading the VCF file for filtering out all sites with no genotyped samples #
#                                                                            #
##############################################################################

# Relative input string for the pyVCF module
relVcfPathString = sys.argv[1]
# Turning it into absolute path
vcfPathString = pathlib.Path(relVcfPathString)
absVcfPath = vcfPathString.resolve(strict=True)
absVcfPath = str(absVcfPath)
# Reading it as VCF file
vcfFile = vcf.Reader(filename=absVcfPath)

#############################################################
#                                                           #
# Creating output file for missing sample filtered variants #
#                                                           #
############################################################

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


# Filtering the SNPs according to the minimum number of reference and alternative allele
for rec in vcfFile:
    # Counter for heterozygous genotypes
    HET = 0
    # Loop over each sample
    for sample in rec.samples:
        # Checking for heterozygosity
        if sample.called and len(set(sample['GT'].split('/'))) == 2:
            HET += 1
            break # the first heterozygous genotype fulfills the filtering criterium. Jump directly to the next variant
    if HET == 1:
        vcfOutFile.write_record(rec)