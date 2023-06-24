import logging
from shlex import split as shlexsplit
from subprocess import Popen, PIPE, STDOUT,CalledProcessError

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
