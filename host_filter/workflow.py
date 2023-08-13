import snakemake
import yaml

class workflow:
	def __init__(self):
		self.config_file='config.yaml'
		self.snakefile='host_filter/etc/Snakefile'
		self.cluster_config='cluster_config.json'
	
	def unlock(self):
		"""
		Unlocks snakemake directory in event of failed run leaving locks
		intact

		Required arguments: None

		Returns: None
		"""

		status = snakemake.snakemake(
			self.snakefile, 
			config={'mapq':None, 'kraken_db':None, 'metaphlan_db': None, 'results': None},
			configfiles = [self.config_file],
			unlock=True
		)

		
	def map(self, verbose):

		"""
		Carries out BWA mapping of paired reads located in `fastq` directory

		Required arguments:
			verbose(bool): Generate verbose output

		Returns:
			None
		"""

		with open(self.config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']
		
		status = snakemake.snakemake(
			self.snakefile, 
			targets = ['mappings'], 
			printshellcmds = False,
			printreason = True,
			nodes = 100, 
			latency_wait = 60, 
			config={'mapq':None, 'kraken_db':None, 'metaphlan_db': None, 'results': None},
			configfiles = [self.config_file],
			cluster_config = self.cluster_config,
			drmaa = drmaa, 
			keepgoing = True, 
			jobname = "map_{jobid}",
			shadow_prefix = '/tmp', 
			use_conda = True,
			conda_frontend='mamba',
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

		with open(self.config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']

		status = snakemake.snakemake(
	 	self.snakefile,
			targets = ['mapq_filter'],
			printshellcmds = False,
			printreason = True,
			nodes = 100,
			latency_wait = 60,
			config={'mapq':mapq, 'kraken_db':None, 'metaphlan_db':None, 'results': None},
			configfiles = [self.config_file],
			cluster_config = self.cluster_config,
			drmaa = drmaa,
			keepgoing = True,
			jobname = "filter_mapq_{wildcards.sample}_{wildcards.mapq}_{jobid}",
			shadow_prefix = '/tmp',
			verbose=verbose,
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
		with open(self.config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']

		mapq=list(range(start_mapq, end_mapq, increment))
		if verbose:
			print(f'Filtering MAPQ {mapq}')

		status = snakemake.snakemake(
			self.snakefile,
			targets = ['mapq_filter'],
			printshellcmds = False,
			printreason = True,
			nodes = 100,
			latency_wait = 60,
			config={'mapq':mapq, 'kraken_db':None, 'metaphlan_db':None, 'results': None},
			configfiles = [self.config_file],
			cluster_config = self.cluster_config,
			drmaa = drmaa,
			keepgoing = True,
			jobname = "filter_mapq_{wildcards.sample}_{wildcards.mapq}_{jobid}",
			shadow_prefix = '/tmp',
			verbose=verbose,
			lint=False
		)
		if not status:
			raise Exception('Snakemake failed filtering MAPQ {mapq}')

	def metaphlan(self, db, verbose):

		"""
		Classifies filtered fastq files using metaphlan

		Required parameters: 
			verbose (bool): verbose reporting
			db(str): database name to search

		Returns:
			None
		"""

		with open(self.config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']

		if verbose:
			print('Classifying with Metaphlan')
			
		status = snakemake.snakemake(
		self.snakefile,
			targets = ['metaphlan'],
			printshellcmds = True,
			printreason = True,
			nodes = 100,
			latency_wait = 60,
			config={'mapq':'filtered', 'kraken_db':None, 'metaphlan_db':db, 'results': None},
			configfiles = [self.config_file],
			cluster_config = self.cluster_config,
			drmaa = drmaa,
			keepgoing = True,
			jobname = "metaphlan_{wildcards.sample}_{wildcards.state}_{jobid}",
			shadow_prefix = '/tmp',
			use_conda=True,
			conda_frontend='mamba',
			verbose=verbose,
			lint=False,
			summary=False,
			dryrun=False
		)
		if not status:
			raise Exception('Snakemake failed running metaphlan}')

	def kraken(self, db, verbose):

		"""
		Classifies filtered fastq files using kraken

		Required parameters: 
			db: (str): database name to search
			verbose (bool): verbose reporting

		Returns:
			None
		"""

		with open(self.config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']

		if verbose:
			print(f'Classifying with Kraken using {db}')
			
		status = snakemake.snakemake(
		self.snakefile,
			targets = ['kraken'],
			printshellcmds = True,
			printreason = True,
			nodes = 100,
			latency_wait = 60,
			config={'mapq':'filtered', 'kraken_db':db, 'metaphlan_db':None, 'results': None},
			configfiles = [self.config_file],
			cluster_config = self.cluster_config,
			drmaa = drmaa,
			keepgoing = True,
			shadow_prefix = '/tmp',
			use_conda=True,
			conda_frontend='mamba',
			verbose=verbose,
			lint=False,
			summary=False,
			dryrun=False
		)
		if not status:
			raise Exception('Snakemake failed running kraken}')

	def analyse(self, results, mapq, verbose):

		"""
		Carries out downstream analysis on obtained results

		Required parameters: 
			results: (str): results set to search
			mapq: (int): Mapq set to analyse
			verbose (bool): verbose reporting

		Returns:
			None
		"""

		with open(self.config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']

		if verbose:
			print(f'Analysing {mapq} results set...')
			
		status = snakemake.snakemake(
		self.snakefile,
			targets = ['analyse'],
			printshellcmds = True,
			printreason = True,
			nodes = 100,
			latency_wait = 60,
			config={'mapq':mapq, 'kraken_db':results},
			configfiles = [self.config_file],
			cluster_config = self.cluster_config,
			drmaa = drmaa,
			keepgoing = True,
			jobname = "analyse_{wildvards.results}_{jobid}",
			shadow_prefix = '/tmp',
			use_conda=False,
			conda_frontend='mamba',
			verbose=verbose,
			lint=False,
			summary=False,
			dryrun=False
		)
		if not status:
			raise Exception('Snakemake failed running analyse}')