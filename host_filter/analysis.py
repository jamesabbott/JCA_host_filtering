from glob import glob
from host_filter import common,genome
from biom import load_table
from biom.table import Table
from pathlib import Path
from tqdm import tqdm
import os
import pandas as pd
import re

class analysis:
	def __init__(self):
		self.config=common.read_config()
		self.db_obj=genome.database()
		pass

	def merge_biom(self, results, mapq, verbose):
		"""
		Merges all available biom files for a mapping quality threshold

		Required parameters:
			results(str): results set
			mapq(str): mapping quality threshold
			verbose(bool): Report verbose output

		Returns:
			merged(str): path to merged biom file	
		"""
		if not os.path.exists(f'analysis/{results}/{mapq}.biom'):
			bioms=glob(f"kraken/{results}/{mapq}/*biom")	
			uber_table=None

			with tqdm(total=len(bioms)) as pbar:
				pbar.set_description(f'Merging {mapq} biom files...')
				for biom in bioms:
					table=load_table(biom)
					if not uber_table:
						uber_table=table
					else:
						uber_table=uber_table.merge(table)

					pbar.update()

			with open(f'analysis/{results}/{mapq}.biom','w') as out_fh:
				out_fh.write(uber_table.to_json(generated_by='host_filter.py'))

		return(f'analysis/{results}/{mapq}.biom')
	
	# def combine_kraken_results(self, results, mapq):
	# 	"""
	# 	Generates summary text file from kraken results summarising
	# 	per-sample host and kingdom representation in each sample.

	# 	Required params: 
	# 	results(str): results set

	# 	Returns:
	# 	None
	# 	"""
	# 	db_info=self.db_obj.get_db_info(self.config['database'])
	# 	species=db_info.get('species')
	# 	superkingdom=db_info.get('superkingdom')

	# 	res_files=glob(f'kraken/{results}/{mapq}/*report.txt')
	# 	results=list()
	# 	print(species)
	# 	for file in res_files:
	# 		kraken_df=pd.read_csv(file, sep="\t",header=None, 
	# 			names=('frag_percent','frag_clade','frag_direct','rank','taxid','sci_name'))
	# 		print(kraken_df[kraken_df['sci_name'].str contains(species)])
			
	# 		import sys
	# 		sys.exit(0)
			
	def match_re(self,sample_data, regex,string):
		match=regex.match(string)
		if match:
			sample_data[match.group(2)]=match.group(1)

	def summarise_mapping_results(self):
		"""
		For a given mapq cutoff, generates per-sample summary tables
		of mapping statistics from samtools flagstat outputs

		Required parameters:
			None

		Returns:
			summary(pd.DataFrame): summary table of values
		"""
		total_re=re.compile(r'(^[\d]+) \+ [\d]+ in (total)')
		mapped_re=re.compile(r'(^[\d]+) \+ [\d]+ (mapped)')
		primary_re=re.compile(r'(^[\d]+) \+ [\d]+ (primary) mapped')
		paired_re=re.compile(r'(^[\d]+) \+ [\d]+ properly (paired)')

		results=dict()
		files=glob('mappings/*flagstat')
		for file in files:
			sample=Path(file).name.removesuffix('.flagstat')
			sample_data=dict()

			with open(file,'r') as fh:
				lines=fh.readlines()
				for line in lines:
					self.match_re(sample_data, mapped_re, line)
					self.match_re(sample_data, total_re, line)
					self.match_re(sample_data, primary_re, line)
					self.match_re(sample_data, paired_re, line)

			results[sample]=sample_data

		mapping_summary=pd.DataFrame.from_dict(results, orient='index')
		mapping_summary.to_csv(f'analysis/mapping_summary.txt',sep="\t", header=True, index=True)

	def analyse(self, results, mapq, verbose):
		"""
		Carries out downstream analysis on a set of analyse results

		Initially merges biom files, then calls R script to carry
		out comparative analysis

		Required paramaters: 
			results(str): Name of results set
			mapq(str): mapq data to process
			verbose(bool): Display verbose output

		Returns:
			None
		"""
		method=results.split('_')[0]
		os.makedirs(f'analysis/{results}',exist_ok=True)

		bioms=list()
		merged=self.merge_biom(results, mapq, verbose)
		bioms.append(merged)
		
		mapping_summary=self.summarise_mapping_results()

		#if method == 'kraken':
		#self.combine_kraken_results(results, mapq)

		return()
