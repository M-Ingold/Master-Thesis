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


#####################################################################################
#                                                                                   #
# Reading the mean read depth file for filtering for non-called sites and max depth #
#                                                                                   #
#####################################################################################

# Relative input string for the pandas module
relDepPathString = "out.ldepth.mean"
# Turning it into absolute path
DepPathString = pathlib.Path(relDepPathString)
absDepPath = DepPathString.resolve(strict=True)
absDepPath = str(absDepPath)
# Reading it as VCF file
meanDepthFrame = pd.read_csv(absDepPath, sep="\t")


#########################################################
#                                                       #
# Creating output file for max. depth filtered variants #
#                                                       #
#########################################################

# Relative input string for the pyVCF module
#relVcfOutPathString = "../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth.vcf"
relVcfOutPathString = sys.argv[2]
# Turning it into absolute path
vcfOutPathString = pathlib.Path(relVcfOutPathString)
absVcfOutPath = vcfOutPathString.resolve(strict=True)
absVcfOutPath = str(absVcfOutPath)
# Creating the output object using the output string
vcfOutFile = vcf.Writer(open(absVcfOutPath, 'w'), vcfFile)


#####################################################################################
#                                                                                   #
# Actual filtering and writing variants that passed to the  #
#                                                                                   #
#####################################################################################


# Setting the 0.95 quantile as threshold for maximum mean depth
maxThreshold = np.quantile(meanDepthFrame["MEAN_DEPTH"], 0.95) # ~769.6
# Iterating over each entry in the depth data frame (each entry corresponds to one variant)
for i in range(meanDepthFrame.shape[0]):
    rec = next(vcfFile)
    # Checking the concordance between the variant's position and the position in the depth data frame
    if rec.CHROM == meanDepthFrame['CHROM'][i] and rec.POS == meanDepthFrame['POS'][i]:
        # Writing selected variant to the output file, if its mean sequencing depth is below the upper threshold
        if meanDepthFrame["MEAN_DEPTH"][i] <= maxThreshold:
            vcfOutFile.write_record(rec)