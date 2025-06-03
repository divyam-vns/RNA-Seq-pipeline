#!/bin/bash -l

module load STAR/2.7.10a
module load samtools/1.15.1

## STAR script
STAR=STAR
BASEDIR=/scratch16/lflorea1/projects/Coursera/M1/
DATADIR=${BASEDIR}/Data/FASTQ
WORKDIR=${BASEDIR}/Star
STAR2IDX=/scratch4/lflorea1/shared/genomes/mm39/star


for i in SRR3734796 SRR3734797 SRR3734798
do
  time \
  $STAR --limitBAMsortRAM 50000000000 --runThreadN 6 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq.gz --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix ${WORKDIR}/${i}_ &> ${WORKDIR}/${i}.star.log &
done

wait 

for i in SRR3734820 SRR3734824 SRR3734825
do
  time \
  $STAR --limitBAMsortRAM 50000000000 --runThreadN 6 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq.gz --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix ${WORKDIR}/${i}_ &> ${WORKDIR}/${i}.star.log &
done


wait

exit;

