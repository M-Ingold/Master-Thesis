#!/bin/bash

# input a target VCF to be diploidized, converted to BED format and analysed by STRUCTURE, 
# using as many clusters as specified. In this case 10

IN_FILE=$1
NAME=$(basename $1 .vcf)
FOLDER=../../data/diploid_VCF # diploid VCF is saved in a separate folder

# target file gets diploidized
sh ../diploid_variant_call/tetraploid2diploid.sh $IN_FILE

# target file is converted to BED format
NEW_FILE=$FOLDER/${NAME}_diploidized.vcf
~/genetools/plink_linux_x86_64_20210606/plink --vcf $NEW_FILE --recode --make-bed --out $NEW_FILE

BED_FILE=$NEW_FILE.bed

for K in 1 2 3 4 5 6 7 8 9 10; \
do ~/genetools/admixture_linux-1.3.0/admixture --cv $BED_FILE $K | tee log${K}.out; done

# display lowest error K for each number of K
echo "Pick lowest error:"
x=$(grep -h CV log*.out)
echo $x