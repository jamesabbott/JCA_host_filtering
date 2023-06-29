import pysam
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
		outbam=f'{sample}.mapq_{mapq}_filtered.bam'

		makedirs(f'filtered_fastq/mapq_{mapq}',exist_ok=True)
		makedirs(f'filtered_bams/mapq_{mapq}',exist_ok=True)
			
		sam=pysam.AlignmentFile(Path(f'mappings/{bam}'),'rb')
		outsam=pysam.AlignmentFile(Path(f'filtered_bams/mapq_{mapq}/{outbam}'), "wb", template=sam)

		for read in sam.fetch(until_eof=True):
			if read.is_unmapped:
				outsam.write(read)
			elif read.is_paired and read.mapping_quality < mapq:
				outsam.write(read)

		sam.close
		outsam.close

		pysam.index('-c',f'filtered_bams/mapq_{mapq}/{outbam}')
		pysam.fastq('-N', '-1', f'filtered_fastq/mapq_{mapq}/{sample}_unmapped_1.fq.gz',
						  '-2', f'filtered_fastq/mapq_{mapq}/{sample}_unmapped_2.fq.gz',
						  '-0', '/dev/null', f'filtered_bams/mapq_{mapq}/{outbam}')

