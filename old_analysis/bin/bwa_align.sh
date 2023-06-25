#!/bin/env bash

#$ -pe smp 24
#$ -mods l_hard h_vmem 32G
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -cwd 
#$ -t 1-12:1

readarray -t samples < <(ls -1 fastq|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]}

mkdir -p bwa_alignments

cp -rv bwa/* $TMPDIR/
cp -v fastq/${sample}* $TMPDIR
bwa mem -t 24 -M $TMPDIR/MorexV3 $TMPDIR/${sample}_1.fq.gz $TMPDIR/${sample}_2.fq.gz |\
	samtools view -b -o $TMPDIR/${sample}.bam
samtools sort -@ 24 --write-index -o $TMPDIR/${sample}.sorted.bam $TMPDIR/${sample}.bam

cp -v ${TMPDIR}/${sample}.sorted* bwa_alignments
