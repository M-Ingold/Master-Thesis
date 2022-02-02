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

IN_DEPTH=$OUT_SAMPLES
OUT_DEPTH=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth.vcf

IN_BLANK=$OUT_DEPTH
OUT_BLANK=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked.vcf

IN_BIAS=$OUT_BLANK
OUT_BIAS=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read.vcf

IN_HET=$OUT_BIAS
OUT_HET=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het.vcf

OUT_SNPs=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_SNPs.vcf
OUT_MNPs=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_MNPs.vcf
OUT_INDELs=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_INDELs.vcf
OUT_COMPLEX=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_COMPLEX.vcf

IN_BIALLEL=$OUT_SNPs
OUT_BIALLEL=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs.vcf

IN_SAMPLE_DEPTH=$OUT_BIALLEL
OUT_SAMPLE_DEPTH=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_dp.vcf

IN_MISSING=$OUT_BIALLEL
OUT_MISSING=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked.vcf

IN_MAF=$OUT_MISSING
OUT_MAF=../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs_blanked_MAF.vcf

#QUAL >= 30
#sh filter_VCF_GWAS_QUAL.sh $IN_QUAL $OUT_QUAL

#filter out double genotypes: P2-AT, P2-LI, P1-AG, P01-BEL, P2-F12-172-DE, P2-G12-173-SO, P2-H12-179-AG 
#and bad coverage sample BIR-2. List created by hand to-do P2 AG statt P1
#bcftools view --samples-file samples.txt $IN_SAMPLES > $OUT_SAMPLES

# Create statistics for average depth for subsequent filtering of max depth
#IN_PATH_SUBSET=../../data/vcf_GWAS_input
#IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30.vcf
#OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS
#mkdir -p $OUT_STAT_PATH
#OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_135-samples_Biall_SNPs_QUAL-30
# Extracting information from the subsetted SNP VCF file
#vcftools --vcf $OUT_SAMPLES --site-mean-depth
#vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
#vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME

# Depth filter
#touch $OUT_DEPTH
#python3 filter_VCF_GWAS_depth_per_site.py $IN_DEPTH $OUT_DEPTH


# Blanking "bad" genotypes and filtering variants with too few genotyped samples !!!filter criteria commented out since they only worked for biallelic SNPs so far, to do later!!!
#touch $OUT_BLANK
#python3 filter_VCF_GWAS_Missing.py $IN_BLANK $OUT_BLANK


# Min. 1 read per strand per variant
#sh filter_VCF_GWAS_StrandBias.sh $IN_BIAS $OUT_BIAS

# Min. 1 ALT and 1 REF allele in the population at this site
#touch $OUT_HET
#python3 filter_VCF_GWAS_HET.py $IN_HET $OUT_HET

#Separating Variant Types into unique files
#sh separate_VCF_types.sh $OUT_HET $OUT_SNPs $OUT_MNPs $OUT_INDELs $OUT_COMPLEX

#Filtering for biallelic SNPs
#touch $OUT_BIALLEL
#python3 filter_VCF_GWAS_biallelic_SNPs.py $IN_BIALLEL $OUT_BIALLEL

##filter for sample depth?
##gzip ../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs.vcf 
##touch $OUT_SAMPLE_DEPTH
##python3 vcf_filter_from_pyVCF.py ../../data/VCF/freebayes_261_samples_chr01-12_QUAL_30_depth_0.9_blanked_1_read_het_biallelic_SNPs.vcf.gz dps --depth-per-sample 61 > $OUT_SAMPLE_DEPTH

#apply full blanking script
#touch $OUT_MISSING
#python3 filter_VCF_GWAS_Missing_full_script.py $IN_MISSING $OUT_MISSING

 ~/genetools/bcftools-1.14/bcftools view --min-af 0.05 --max-af 0.95 --exclude-uncalled $IN_MAF > $OUT_MAF