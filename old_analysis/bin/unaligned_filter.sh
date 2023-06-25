#!/bin/env bash

# filter to extract unaligned reads as was done for the NT paper

#$ -pe smp 24
#$ -mods l_hard h_vmem 32G
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -cwd 
#$ -t 1-12:1

readarray -t samples < <(ls -1 fastq|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]}

mkdir -p filtered_fastq/morex3_unaligned

cp -rv bwa_alignments/${sample}* $TMPDIR/
samtools fastq -@ 24 -f 12 -N -1 $TMPDIR/${sample}_unmapped_1.fq.gz -2 $TMPDIR/${sample}_unmapped_2.fq.gz \
	-0 /dev/null $TMPDIR/${sample}.sorted.bam

cp -v ${TMPDIR}/${sample}_unmapped*.gz filtered_fastq/morex3_unaligned