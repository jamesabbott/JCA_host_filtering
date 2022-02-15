#!/bin/env bash

#$ -pe smp 24
#$ -mods l_hard mfree 128G
#$ -mods l_hard local_free 200G
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -cwd
#$ -l h=c6320*

set -e

# Runs kraken2 against paired fastq reads, where read filenames start with the sample followed by a '_'.

# This script is designed to run a single cluster job - the overhead of copying the database to the node vs the 
# kraken runtime makes this more effecient than running a separate task per-sample
# The database is copied to the $TMPDIR of the execution host 
# which has local disk storage available on execution nodesaccessed via $TMPDIR, 

usage="Usage: $0 -i /path/to/input/fastq/directory -o /path/to/output/directory -d path/to/kraken_db" 

while getopts "i:o:d:" opt; do
  case "${opt}" in
	i )
		in_dir=${OPTARG}
		;;
	o )
	  out_dir=${OPTARG}
	  ;;
	d )
		db=${OPTARG}
	  ;;
	\? )
		echo ${usage}
	  ;;
	: )
		echo "Invalid option: $OPTARG requires an argument" 1>&2
	;;
  esac
done

if [[ ! -e "${in_dir}" ]]; then
	echo "Specificed input directory (${indir}) does not exist"
	exit
fi

if [[ ! -e "${db}" ]]; then
	echo "Specificed kraken database (${db}) does not exist"
	exit
fi

mkdir -p ${out_dir}

readarray -t samples < <(ls $in_dir|cut -f1 -d_|uniq)


df -h

mkdir $TMPDIR/kraken_db
cp -v ${db}/*k2d $TMPDIR/kraken_db
cp -v ${in_dir}/* $TMPDIR
ls -l $TMPDIR
ls -l $TMPDIR/kraken_db

for sample in ${samples[@]}; do
	echo sample=${sample}
	readarray -t files < <(ls -1 $in_dir/${sample}*|xargs -i basename {})
	sample=$(echo $sample|sed 's/.filtered//')
	kraken2 --db $TMPDIR/kraken_db --threads 24 --use-names --report ${out_dir}/${sample}.report.txt \
		--output ${out_dir}/${sample}.txt --gzip-compressed --paired ${TMPDIR}/${files[0]} ${TMPDIR}/${files[1]}
	#bracken -d ${TMPDIR}/kraken_db -i ${out_dir}/${sample}.report.txt -o ${out_dir}/${sample}.bracken -r 150 
done
