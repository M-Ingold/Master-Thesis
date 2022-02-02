#!/bin/bash

################################################################
#                                                              #
# Master script for applying all selected filters sequencially #
#                                                              #
################################################################

# In and output for QUAL filtering
IN_QUAL=../../data/VCF/freebayes_269_samples_chr01-12.vcf
OUT_QUAL=../../data/VCF/freebayes_269_samples_chr01-12_QUAL_30.vcf

IN_SAMPLES=$OUT_QUAL
OUT_SAMPLES=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30.vcf

OUT_SNPs=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs.vcf
OUT_MNPs=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_MNPs.vcf
OUT_INDELs=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_INDELs.vcf
OUT_COMPLEX=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_COMPLEX.vcf

IN_BIAS=$OUT_SNPs
OUT_BIAS=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs_1_read.vcf

IN_HET=$OUT_BIAS
OUT_HET=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_SNPs_1_read_het.vcf

IN_BIALLEL=$OUT_HET
OUT_BIALLEL=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs.vcf

IN_MISSING=$OUT_BIALLEL
OUT_MISSING=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked.vcf

IN_DEPTH=$OUT_MISSING
OUT_DEPTH=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth.vcf

IN_MAF=$OUT_DEPTH
OUT_MAF=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_1_read_het_biallelic_SNPs_blanked_depth_MAF.vcf

#QUAL >= 30
#sh filter_VCF_GWAS_QUAL.sh $IN_QUAL $OUT_QUAL

#filter out double genotypes: P2-AT, P2-LI, P1-AG, P01-BEL, P2-F12-172-DE, P2-G12-173-SO, P2-H12-179-AG 
#and bad coverage sample BIR-2. List created by hand
bcftools view --samples-file samples.txt $IN_SAMPLES > $OUT_SAMPLES

#Separating Variant Types into unique files
sh separate_VCF_types.sh $OUT_SAMPLES $OUT_SNPs $OUT_MNPs $OUT_INDELs $OUT_COMPLEX

# Min. 1 read per strand per variant
touch $OUT_BIAS
sh filter_VCF_GWAS_StrandBias.sh $IN_BIAS $OUT_BIAS

# Min. 1 ALT and 1 REF alleles or 2 ALT and 0 REF alleles in the population at this site
touch $OUT_HET
python3 filter_VCF_GWAS_HET.py $IN_HET $OUT_HET

#Filtering for biallelic SNPs
touch $OUT_BIALLEL
python3 filter_VCF_GWAS_biallelic_SNPs.py $IN_BIALLEL $OUT_BIALLEL

# Blanking "bad" genotypes and filtering variants with too few genotyped samples
touch $OUT_MISSING
python3 filter_VCF_GWAS_Missing.py $IN_MISSING $OUT_MISSING

# Create statistics for average depth for subsequent filtering of max depth
#IN_PATH_SUBSET=../../data/vcf_GWAS_input
#IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30.vcf
#OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS
#mkdir -p $OUT_STAT_PATH
#OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_135-samples_Biall_SNPs_QUAL-30
# Extracting information from the subsetted SNP VCF file
vcftools --vcf $IN_DEPTH --site-mean-depth
#vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
#vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME

# Max depth filter
touch $OUT_DEPTH
python3 filter_VCF_GWAS_depth_per_site.py $IN_DEPTH $OUT_DEPTH

bcftools view --min-af 0.05 --max-af 0.95 --exclude-uncalled $IN_MAF > $OUT_MAF