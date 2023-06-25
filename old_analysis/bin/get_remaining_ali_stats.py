#!/bin/env python

#$ -cwd
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -pe smp 8

import pysam
import pandas as pd
import shutil
import os
from multiprocessing import Pool

def get_stats(chr):
#\for read in sam.fetch(until_eof=True):    
	print('Processing {}'.format(chr))
	max_reads=100000
	read_count=0

	alignment_stats=pd.DataFrame(columns=['read','cigar','mapq','length'])
	sam=pysam.AlignmentFile('{}/alignments.bam'.format(os.environ['TMPDIR']),'rb')
	for read in sam.fetch(chr,multiple_iterators=True):
		if read.query_name in barley_read_ids.values:
			stats={
				'read':read.query_name,
				'cigar':read.get_cigar_stats(),
				'mapq':read.mapping_quality,
				'length':read.query_alignment_length
			}
			aliignment_stats=alignment_stats.append(stats,ignore_index=True)
			read_count+=1

		if read_count==max_reads:
			print('Writing {} results...'.format(chr))
			alignment_stats.to_csv('mapq30_remaining_barley_stats_{}.txt'.format(chr),sep="\t",index=False)
	sam.close()

kraken_report=pd.read_csv('kraken/nt_mapq30_minlen100/2000.txt',sep="\t",header=None)
kraken_report=kraken_report.loc[kraken_report[2]=='Hordeum vulgare (taxid 4513)']
barley_read_ids=kraken_report[1]

shutil.copyfile('bwa_alignments/2000.mapq_30_minlen_100_filtered.fixed.bam','{}/alignments.bam'.format(os.environ['TMPDIR']))
shutil.copyfile('bwa_alignments/2000.mapq_30_minlen_100_filtered.fixed.bam.csi','{}/alignments.bam.csi'.format(os.environ['TMPDIR']))
sam=pysam.AlignmentFile('{}/alignments.bam'.format(os.environ['TMPDIR']),'rb')
chrs=sam.header.references
print(chrs)
sam.close()

get_stats('chr1H')
#pool=Pool(processes=8)
#pool.map(get_stats,chrs)
#pool.close
#pool.join


