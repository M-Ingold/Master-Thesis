#!/bin/bash

IN_FILE=$1
NAME=$(basename $1 .vcf)
FOLDER=../../data/diploid_VCF
OUT_FILE=${NAME}_diploidized.vcf

cp $IN_FILE $FOLDER

sed -i 's+1/1/1/1+1/1+g' $FOLDER/$NAME.vcf
sed -i 's+0/1/1/1+0/1+g' $FOLDER/$NAME.vcf
sed -i 's+0/0/1/1+0/1+g' $FOLDER/$NAME.vcf
sed -i 's+0/0/0/1+0/1+g' $FOLDER/$NAME.vcf
sed -i 's+0/0/0/0+0/0+g' $FOLDER/$NAME.vcf

mv $FOLDER/$NAME.vcf $FOLDER/$OUT_FILE