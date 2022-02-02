#!/bin/bash

BAM_FOLDER=../../data/alignment

# this only outputs summary stats, not per base coverage
##~/genetools/samtools-1.14/samtools coverage --min-read-len 1 --depth 0 --output $BAM_FOLDER/samtools_coverage.txt $BAM_FOLDER/all_samples.bam 

bedtools genomecov -dz -ibam $BAM_FOLDER/all_samples.bam -g ../../References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa > $BAM_FOLDER/bedtools_coverage.txt

# subset by chromosome
i=1

while [ $i -lt 10 ]; do
	grep chr0$i $BAM_FOLDER/bedtools_coverage.txt > $BAM_FOLDER/coverage_chr0$i.txt
	i=$(($i+1))
done

while [ $i -lt 13 ]; do
	grep chr$i $BAM_FOLDER/bedtools_coverage.txt > $BAM_FOLDER/coverage_chr$i.txt
	i=$(($i+1))
done