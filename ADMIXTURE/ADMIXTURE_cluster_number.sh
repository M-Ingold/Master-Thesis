#!/bin/bash
# PED File not working for some reason

BED_FILE=../../data/diploid_VCF/diploid.bed

for K in 1 2 3 4 5 6 7 8 9 10; \
do ~/genetools/admixture_linux-1.3.0/admixture --cv $BED_FILE $K -jN 11| tee log${K}.out; done

# pick lowest error K
echo "Pick lowest error:"
x=$(grep -h CV log*.out)
echo $x