PATH_IN_BAMS=../data/alignment
PATH_BAM_LIST=bam_list_for_bamaddrg.txt

# Creating a space separated list of input BAM files of the form "-b $BAM_FILE_PATH -s $SAMPLE_NAME
for folder in ${PATH_IN_BAMS}/Sample*; do
	sample=$(basename $folder)
	bamFile=$(realpath ${folder}/${sample}_sorted.bam)
	echo -n "-b ${bamFile} -s ${sample} " >> $PATH_BAM_LIST
done


TOTAL_SAMPLES=$(ls -l $PATH_IN_BAMS | grep -c ^d)

BIG_BAM_FILE=$PATH_IN_BAMS/all_samples.bam

BAMADDRG_INPUT=$(head -n 1 ${PATH_BAM_LIST}) # This is necessary for storing the file content in a variable, as expected from bamadrrg

~/genetools/bamaddrg/bamaddrg --clear $BAMADDRG_INPUT > $BIG_BAM_FILE

samtools index $BIG_BAM_FILE

rm $PATH_BAM_LIST