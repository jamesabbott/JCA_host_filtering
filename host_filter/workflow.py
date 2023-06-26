import snakemake
import yaml

class workflow:
	def __init__(self):
		pass
		
	def run(self, action):
		print(f'running {action}')

		config_file='config.yaml'
		snakefile='host_filter/etc/Snakefile'
		cluster_config='cluster_config.json'

		with open(config_file,'r') as in_fh:
			config=yaml.safe_load(in_fh)

		drmaa=config['drmaa']
		
		status = snakemake.snakemake(
			snakefile, 
			targets = [action], 
			printshellcmds = False,
			printreason = True,
			nodes = 100, 
			latency_wait = 60, 
			configfiles = [config_file],
			cluster_config = cluster_config,
			drmaa = drmaa, 
			keepgoing = True, 
			jobname = "host_filter_{jobid}",
			shadow_prefix = '/tmp', 
			use_conda = True
		)

		if status: 
			return 0
		return 1