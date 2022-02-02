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
#############################################################

# Relative input string for the pyVCF module
#relVcfOutPathString = "../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth_MissingsSamples-50.vcf"
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

# Filtering only for biallelic SNPs
for rec in vcfFile:
    # Check for SNPs that have only one alternative allele besides the reference allele
    if len(rec.INFO["AO"]) == 1:
        vcfOutFile.write_record(rec)
    # In case of absence of a reference allele in the population, check for SNPs that have exactly two alternative alleles
    elif rec.INFO["RO"] == 0 and len(rec.INFO["AO"]) == 2:
        vcfOutFile.write_record(rec)