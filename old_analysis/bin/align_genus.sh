#!/bin/env bash

#$ -pe smp 24
#$ -mods l_hard h_vmem 32G
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -cwd 
#$ -t 1-22:1

readarray -t samples < <(ls -1 per_genus/fastq|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]}

mkdir -p per_genus/bwa_alignments

cp -rv bwa/* $TMPDIR/
cp -v per_genus/fastq/${sample}* $TMPDIR
bwa mem -t 24 -M $TMPDIR/MorexV3 $TMPDIR/${sample}_1.fq.gz $TMPDIR/${sample}_2.fq.gz |\
	samtools view -b -o $TMPDIR/${sample}.bam
samtools sort -@ 24 --write-index -o $TMPDIR/${sample}.sorted.bam $TMPDIR/${sample}.bam
samtools flagstat $TMPDIR/${sample}.sorted.bam > $TMPDIR/${sample}.sorted.flagstat

cp -v ${TMPDIR}/${sample}.sorted* per_genus/bwa_alignments