from host_filter import common
from json import load
from os import makedirs,environ
from pathlib import Path
from tabulate import tabulate
from yaspin import yaspin
import tarfile
import requests

"""
Methods relating to classification databases, execution etc.
"""

class classify:
	def __init__(self):
		self.db_path=Path(environ['CONDA_PREFIX']) / Path('share/classification_databases')
		self.dbs=self.read_dbs()
	
	def read_dbs(self):
		"""
		Obtains info on configured databases from json file

		Required params: 
			None

		Returns:
			dbs(dict): keyed on database name
		"""
		path=Path.resolve(Path('.'))
		with open(path / Path('host_filter/etc/classify_dbs.json')) as fh:
			dbs=load(fh)
			return(dbs)
	
	def available_dbs(self):
		"""
		Lists databases available for download

		Required parameters: 
			None

		Returns: 
			None
		"""

		output=[['Database','Downloaded']]
		for db in self.dbs:
			have_db='No'

			done_file=self.db_path / Path(f'{db}/{db}.done')
			if done_file.exists():
				have_db='Yes'
			
			output.append([db,have_db])
			
		print(tabulate(output,headers='firstrow',tablefmt='fancy_grid'))

	def download(self,db):
		"""
		Downloads specified database 

		Required parameterse:
			db(str): database name

		Returns: 
			None
		"""

		db_info=self.dbs[db]
		db_path=self.db_path / Path(db)
		makedirs(db_path ,exist_ok=True)
		done_file=db_path / Path(f"{db}.done")

		with yaspin(text=f"Downloading {db} database...", timer=True):
			command=db_info.get('command')
			file=db_info.get('file')
			if command:
				command=command.replace('__DBPATH__', f"{db_path}")
				common.run_command(f"{command}", f'{db} download')
			elif file: 
				r=requests.get(f"{file}",stream=True)
				file_path=Path(file)
				dl_file= db_path / file_path.name
				
				with open(dl_file,'wb') as fh:
					for chunk in r.iter_content(chunk_size=4096):
						fh.write(chunk)

				tar=tarfile.open(dl_file,'r:gz')
				tar.extractall(path=db_path)
				dl_file.unlink()

		done_file.touch()
