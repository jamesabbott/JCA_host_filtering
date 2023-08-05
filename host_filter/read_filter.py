import pysam
import gzip
from Bio import SeqIO
from pathlib import Path
from os import makedirs

class read_filter:
	def __init__(self):
		pass

	def mapq_filter(self, mapq, sample):
		"""
		filters bam file based upon specified mapq threshold

		Note this differs from conventional mapq filtering in that we want
		to retain reads below the threshold, whereas standard tools retain 
		reads above the threshold

		Required parameters:
			mapq(int): mapping quality threshold
			sample(str): sample name

		returns:
			None
		"""
		print(f'sample={sample}')
		bam=f'{sample}.bam'

		filtered_bam=f'{sample}.mapq_{mapq}_filtered.bam'
		mapped_bam=f'{sample}.mapq_{mapq}_mapped.bam'

		makedirs(f'filtered_fastq/mapq_{mapq}',exist_ok=True)
		makedirs(f'filtered_bams/mapq_{mapq}',exist_ok=True)
			
		sam=pysam.AlignmentFile(Path(f'mappings/{bam}'),'rb')
		filtered_outsam=pysam.AlignmentFile(Path(f'filtered_bams/mapq_{mapq}/{filtered_bam}'), "wb", template=sam)
		mapped_outsam=pysam.AlignmentFile(Path(f'filtered_bams/mapq_{mapq}/{mapped_bam}'), "wb", template=sam)

		for read in sam.fetch(until_eof=True):
			if read.is_unmapped:
				filtered_outsam.write(read)
			elif read.is_paired and read.mapping_quality < mapq:
				filtered_outsam.write(read)
			elif read.is_paired:
				mapped_outsam.write(read)

		sam.close
		filtered_outsam.close

		pysam.index('-c',f'filtered_bams/mapq_{mapq}/{filtered_bam}')
		pysam.fastq('-N', '-1', f'filtered_fastq/mapq_{mapq}/{sample}_unmapped_1.fq.gz',
						  '-2', f'filtered_fastq/mapq_{mapq}/{sample}_unmapped_2.fq.gz',
						  '-0', '/dev/null', f'filtered_bams/mapq_{mapq}/{filtered_bam}')
		
		pysam.index('-c',f'filtered_bams/mapq_{mapq}/{mapped_bam}')
		pysam.fastq('-N', '-1', f'filtered_fastq/mapq_{mapq}/{sample}_mapped_1.fq.gz',
						  '-2', f'filtered_fastq/mapq_{mapq}/{sample}_mapped_2.fq.gz',
						  '-0', '/dev/null', f'filtered_bams/mapq_{mapq}/{mapped_bam}')

		unmapped_count=0
		mapped_count=0

		with gzip.open(f'filtered_fastq/mapq_{mapq}/{sample}_unmapped_1.fq.gz','rt') as in_fh:
			records=SeqIO.parse(in_fh, 'fastq')
			for record in records:
				unmapped_count=unmapped_count+1

		with gzip.open(f'filtered_fastq/mapq_{mapq}/{sample}_mapped_1.fq.gz','rt') as in_fh:
			records=SeqIO.parse(in_fh, 'fastq')
			for record in records:
				mapped_count=mapped_count+1

		print(mapped_count)
		with open(f'filtered_fastq/mapq_{mapq}/{sample}_stats.txt','w') as stats_fh:
			stats_fh.write(f"Mapped\t{mapped_count}\n")
			stats_fh.write(f"Unmapped\t{unmapped_count}\n")

