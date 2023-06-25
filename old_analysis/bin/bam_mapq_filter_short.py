#!/bin/env python

#$ -cwd
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID
#$ -t 1-6

"""
Filters bam files for alignments based on alignment length to remove 
false +ve hits for host filtering
"""

import os
import pysam
import argparse
from shutil import copyfile

def main():

	parser=argparse.ArgumentParser("Extracts unaligned reads, and those with a MAPQ below a defined threshold from bam alignments")
	parser.add_argument('--mapq',type=int,help="MAPQ threshold")
	parser.add_argument('--min_alignment_length',type=int,help='Minimum length of alignment to retain')
	args=parser.parse_args()

	task_id=int(os.environ['SGE_TASK_ID'])
	tmpdir=os.environ['TMPDIR']
	samples=['2000','2001','2002','2023','2024','2025']
	sample=samples[task_id-1]
	bam='{}.sorted.bam'.format(sample)
	outbam='{}.mapq_{}_minlen_{}_filtered.bam'.format(sample,args.mapq,args.min_alignment_length)

	try:
		os.mkdir('filtered_fastq/mapq_{}_minlen_{}'.format(args.mapq,args.min_alignment_length))
	except FileExistsError:
		pass
		
	copyfile('bwa_alignments/{}'.format(bam),'{}/{}'.format(tmpdir,bam))
	copyfile('bwa_alignments/{}.csi'.format(bam),'{}/{}.csi'.format(tmpdir,bam))
	
	print('HOSTNAME: {}'.format(os.environ['HOSTNAME']))
	print('SAMPLE: {}'.format(sample))

	sam=pysam.AlignmentFile('{}/{}'.format(tmpdir,bam),'rb')
	outsam=pysam.AlignmentFile('{}/{}'.format(tmpdir,outbam), "wb", template=sam)

	for read in sam.fetch(until_eof=True):
		if read.is_unmapped:
			outsam.write(read)
		elif read.mapping_quality < args.mapq:
			outsam.write(read)
		elif read.query_alignment_length < args.min_alignment_length:
			outsam.write(read)

	sam.close
	outsam.close

	pysam.index('-c','{}/{}'.format(tmpdir,outbam))
	pysam.fastq('-N', '-1', '{}/{}_unmapped_1.fq.gz'.format(tmpdir,sample) , 
		'-2', '{}/{}_unmapped_2.fq.gz'.format(tmpdir,sample), '-0', '/dev/null', 
		'{}/{}'.format(tmpdir,outbam))

	copyfile('{}/{}'.format(tmpdir,outbam),'bwa_alignments/{}'.format(outbam))
	copyfile('{}/{}.csi'.format(tmpdir,outbam),'bwa_alignments/{}.csi'.format(outbam))
	copyfile('{}/{}_unmapped_1.fq.gz'.format(tmpdir,sample),
		'filtered_fastq/mapq_{}_minlen_{}/{}_unmapped_1.fq.gz'.format(args.mapq,args.min_alignment_length, sample))
	copyfile('{}/{}_unmapped_2.fq.gz'.format(tmpdir,sample),
		'filtered_fastq/mapq_{}_minlen_{}/{}_unmapped_2.fq.gz'.format(args.mapq,args.min_alignment_length, sample))


if __name__ == "__main__":
	main()