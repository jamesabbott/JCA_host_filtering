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
	args=parser.parse_args()

	task_id=int(os.environ['SGE_TASK_ID'])
	tmpdir=os.environ['TMPDIR']
	samples=['2000','2001','2002','2023','2024','2025']
	sample=samples[task_id-1]
	bam='{}.sorted.bam'.format(sample)
	outbam='{}.mapq_{}_filtered.bam'.format(sample,args.mapq)

	try:
		os.mkdir('filtered_fastq/mapq_{}_nomq0'.format(args.mapq))
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
		elif read.mapping_quality < args.mapq and read.mapping_quality != 0:
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
		'filtered_fastq/mapq_{}_nomq0/{}_unmapped_1.fq.gz'.format(args.mapq,sample))
	copyfile('{}/{}_unmapped_2.fq.gz'.format(tmpdir,sample),
		'filtered_fastq/mapq_{}_nomq0/{}_unmapped_2.fq.gz'.format(args.mapq,sample))


if __name__ == "__main__":
	main()