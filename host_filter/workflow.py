import snakemake
import yaml

class workflow:
	def __init__(self):
		pass
		
	def map(self, verbose):

		"""
		Carries out BWA mapping of paired reads located in `fastq` directory

		Required arguments:
			verbose(bool): Generate verbose output

		Returns:
			None
		"""

		config_file='config.yaml'
		snakefile='host_filter/etc/Snakefile'
		cluster_config='cluster_config.json'

		with open(config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']
		
		status = snakemake.snakemake(
			snakefile, 
			targets = ['mappings'], 
			printshellcmds = False,
			printreason = True,
			nodes = 100, 
			latency_wait = 60, 
			config = {'mapq': None},
			configfiles = [config_file],
			cluster_config = cluster_config,
			drmaa = drmaa, 
			keepgoing = True, 
			jobname = "map_{jobid}",
			shadow_prefix = '/tmp', 
			use_conda = True,
			verbose=verbose
		)

		if status: 
			return 0
		return 1
	
	def filter(self, mapq, verbose):

		"""
		Discards reads mapping to host genome with a mapping quality greater than specified value

		Required parameters:
			mapq(str): Mapping quality threshold
			verbose(bool): Generate verbose output

		Returns:
			None
		"""

		config_file='config.yaml'
		snakefile='host_filter/etc/Snakefile'
		cluster_config='cluster_config.json'

		with open(config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']
		
		status = snakemake.snakemake(
			snakefile,
			targets = ['mapq_filter'],
			printshellcmds = False,
			printreason = True,
			nodes = 100,
			latency_wait = 60,
			config={'mapq':mapq},
			configfiles = [config_file],
			cluster_config = cluster_config,
			drmaa = drmaa,
			keepgoing = True,
			jobname = "filter_mapq_{mapq}_{jobid}",
			shadow_prefix = '/tmp',
			verbose=verbose
		)

		if status: 
			return 0
		return 1
	
	def range_filter(self, start_mapq, end_mapq, increment, verbose):

		"""
		Discards reads mapping to the host genome using a range of different
		mapping quality thresholds 
		
		Required parameters: 
			start_mapq(int): start of mapping quality range
			end_mapq(int): end of mapping quality range
			increment(int): increment size between filtering runs

		Returns:
			None
		"""

		config_file='config.yaml'
		snakefile='host_filter/etc/Snakefile'
		cluster_config='cluster_config.json'

		with open(config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']

		for mapq in range(start_mapq, end_mapq, increment):
			if verbose:
				print(f'Filtering MAPQ {mapq}')
			
			status = snakemake.snakemake(
			snakefile,
				targets = ['mapq_filter'],
				printshellcmds = False,
				printreason = True,
				nodes = 100,
				latency_wait = 60,
				config={'mapq':mapq},
				configfiles = [config_file],
				cluster_config = cluster_config,
				drmaa = drmaa,
				keepgoing = True,
				jobname = "filter_mapq_{mapq}_{jobid}",
				shadow_prefix = '/tmp',
				verbose=verbose
			)
			if not status:
				raise Exception('Snakemake failed filtering MAPQ {mapq}')