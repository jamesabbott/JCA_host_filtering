import subprocess

def run_command(cmd):
	"""
	Runs a shell command...

	Required params:
		cmd(list): command and it's arguments to be run

	Returns:
		None
	"""

	results=subprocess.run(cmd,check=True,capture_output=True,text=True)
	results.check_returncode()