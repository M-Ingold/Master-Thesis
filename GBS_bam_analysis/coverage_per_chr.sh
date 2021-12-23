#!/bin/bash

BAM_FOLDER=../../data/alignment

#samtools depth $FOLDER/all_samples.bam > samtools_coverage.txt

# subset by chromosome
i=1

while [ $i -lt 10 ]; do
	grep chr0$i $BAM_FOLDER/samtools_coverage.txt > $BAM_FOLDER/coverage_chr0$i.txt
	i=$(($i+1))
done

while [ $i -lt 13 ]; do
	grep chr$i $BAM_FOLDER/samtools_coverage.txt > $BAM_FOLDER/coverage_chr$i.txt
	i=$(($i+1))
done