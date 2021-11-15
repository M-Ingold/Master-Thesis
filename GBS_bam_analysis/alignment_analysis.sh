#!/bin/bash

# interactive mode needed!
# analyse aligned files, input parent folder of folders with bam files

#mkdir ../../analysis/$(basename $BAMFOLDER)/coverage

BAMFOLDER=../../data/alignment/

# For allSamples File
# by default the genome is split into 400 windows in which coverage etc is averaged

qualimap bamqc -bam $BAMFOLDER/all_samples.bam -outdir ../../analysis/$(basename $BAMFOLDER) \
--java-mem-size=50G -gff ../../References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 \
--output-genome-coverage ../../analysis/coverage --paint-chromosome-limits
#687 WINDOWS USED!

#to-do: by sequencing group
#http://qualimap.conesalab.org/doc_html/command_line.html#command-line