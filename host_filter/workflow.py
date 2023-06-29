import snakemake
import yaml

class workflow:
	def __init__(self):
		pass
		
	def map(self):

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
			configfiles = [config_file],
			cluster_config = cluster_config,
			drmaa = drmaa, 
			keepgoing = True, 
			jobname = "map_{jobid}",
			shadow_prefix = '/tmp', 
			use_conda = True
		)

		if status: 
			return 0
		return 1
	
	def filter(self, mapq):

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
			shadow_prefix = '/tmp'
		)

		if status: 
			return 0
		return 1