#!/bin/env bash

# filter bam file by mapping quality

#$ -pe smp 24
#$ -mods l_hard h_vmem 32G
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -cwd 
#$ -t 1-6:1

while getopts "m:" opt; do
  case "${opt}" in
    m )
        MAPQ=${OPTARG}
        ;;
    \? )
        echo ${usage}
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "${MAPQ}" ]]; then
	echo "Usage: $0 -m mapq"
	exit 1
fi

samples=(2000 2001 2002 2023 2024 2025)
sample=${samples[$SGE_TASK_ID-1]}

cp -rv bwa_alignments/${sample}.sorted.* $TMPDIR/
samtools view -@ 24 -q ${MAPQ} -h -b -o $TMPDIR/${sample}.mapq_${MAPQ}_filtered.bam \
	--write-index ${TMPDIR}/${sample}.sorted.bam

cp -v ${TMPDIR}/${sample}.mapq* bwa_alignments
