import logging
from yaml import safe_load
from shlex import split as shlexsplit
from subprocess import Popen, PIPE, STDOUT, CalledProcessError

def read_config():
	"""
	Reads config.yaml file and returns results

	Required parameters:
		None

	Returns:
		config(dict): parsed config data
	"""
	with open('config.yaml','r') as in_fh:
		config=safe_load(in_fh)

	return(config)

def run_command(cmd, name):
	"""
	Runs a shell command, writing stderr and stdout to log files

	Required params:
		cmd(str): command and it's arguments to be run
		name(str): name of task for log files

	Returns:
		None
	"""
	process=Popen(shlexsplit(cmd),stdout=PIPE, stderr=STDOUT)
	with process.stdout:
		for line in iter(process.stdout.readline, b''): 
			logging.info(f'{name}: {line.decode("UTF-8").strip()}')
	return_code = process.wait()
	if return_code:
		raise CalledProcessError(return_code, cmd)
