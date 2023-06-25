#!/bin/env bash

# filter to try to avoid artefactual mappings stripping out non-reference reads

# We'll filter twice:
# * unmapped read pairs (both reads of pair unmapped) (-f 12)
# * mapped reads which are not primary alignments -F 268 (not unmapped F+r - 12), plus ('not not' primary (256)...)

# These then need to be merged, and just in case we have picked up a read-pair twice, dedup them by readid.

#$ -pe smp 24
#$ -mods l_hard h_vmem 32G
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -cwd 
#$ -t 1-12:1

readarray -t samples < <(ls -1 fastq|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]}

mkdir -p filtered_fastq/non_primary

cp -rv bwa_alignments/${sample}* $TMPDIR/
samtools fastq -@ 24 -f 12 -N -1 $TMPDIR/${sample}.unmapped_1.fq.gz -2 $TMPDIR/${sample}.unmapped_2.fq.gz \
	-0 /dev/null $TMPDIR/${sample}.sorted.bam
samtools fastq -@ 24 -F 12 -f 256 -N -1 $TMPDIR/${sample}.mapped_not_primary_1.fq.gz -2 $TMPDIR/${sample}.mapped_not_primary_2.fq.gz \
	-0 /dev/null $TMPDIR/${sample}.sorted.bam

cat ${TMPDIR}/*1.fq.gz >> $TMPDIR/${sample}_all_1.fq.gz
cat ${TMPDIR}/*2.fq.gz >> $TMPDIR/${sample}_all_2.fq.gz

seqkit rmdup -n $TMPDIR/${sample}_all_1.fq.gz -o ${TMPDIR}/${sample}_clean_1.fq.gz
seqkit rmdup -n $TMPDIR/${sample}_all_2.fq.gz -o ${TMPDIR}/${sample}_clean_2.fq.gz

cp -v ${TMPDIR}/${sample}_clean*.gz filtered_fastq/non_primary
