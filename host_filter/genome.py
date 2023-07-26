from glob import glob
from host_filter import common
from json import load
from os import makedirs,remove
from pathlib import Path
from tabulate import tabulate
from yaspin import yaspin
import gzip
import logging
import requests
import shutil

class database:
	def __init__(self):
		self.config=common.read_config()
		self.db_path=self.config['database_path']
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
		with open(path / Path('host_filter/etc/genomes.json')) as fh:
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

		output=[['Database','Downloaded','BWA Indexed']]
		for db in self.dbs:
			db_info=self.dbs[db]
			count=0
			found=0
			bwa_index='No'
			have_db='No'

			for file in db_info['files']:
				count+=1
				if (self.db_path / Path(file['filename']+'.download_complete')).exists():
					found+=1
				if (self.db_path / Path(file['filename']+'.sa')).exists():
					bwa_index='Yes'

			if found==count:
				have_db='Yes'
			
			output.append([db,have_db,bwa_index])
			
		print(tabulate(output,headers='firstrow',tablefmt='fancy_grid'))

	def download(self,db):
		"""
		Downloads files for specified database 

		Required parameterse:
			db(str): database name

		Returns: 
			None
		"""
		db_info=self.dbs.get(db)

		makedirs(self.db_path,exist_ok=True)
		with yaspin(text=f"Downloading {db} database...", timer=True):
			for file in db_info['files']:

				dl_file= self.db_path / Path(file['filename']+'.dl')
				final_file= self.db_path / Path(file['filename'])
				done_file= self.db_path / Path(file['filename']+'.download_complete')

				r=requests.get(file['url'],stream=True)

				with open(dl_file,'wb') as fh:
					for chunk in r.iter_content(chunk_size=4096):
						fh.write(chunk)

				if file['compressed']=='True':
					with gzip.open(dl_file,'rb') as in_fh:
						with open(final_file,'wb') as out_fh:
							shutil.copyfileobj(in_fh, out_fh)
						remove(dl_file)
				else:
					shutil.move(dl_file,final_file)
				done_file.touch()

	def bwa_index(self,db):
		"""
		Generates a BWA index for the specified database

		Required params: 
			db(str): database name
		
		Returns:
			None
		"""

		db_info=self.dbs.get(db)

		with yaspin(text=f"BWA indexing {db} database...", timer=True):
			for file in db_info['files']:
				if file['type']=='genome':
					filename=file['filename']
					break

			genome_file=self.db_path / Path(filename)
			cmd=f'bwa index {genome_file}'
			common.run_command(cmd,f'bwa_index_{db}')

	def clean(self,db):
		"""
		Removes files relating to specified database

		Required parameters:
			db(str): database name

		Returns:
			None
		"""
		db_info=self.dbs.get(db)
		
		with yaspin(text=f"Cleaning up {db} database..."):
			logging.info(f'{db}: Database removed...')
			for file in db_info['files']:
				search_path=self.db_path / Path(file['filename']+'*')
				db_files=glob(str(search_path))
				for db_file in db_files:
					remove(db_file)

	def get_db_index_files(self,db):
		"""
		Creates list of BWA index files for a specified database

		Required params: 
			db(str): database name

		Returns:
			db_files(list): list of BWA index files for database
		"""
		db_info=self.dbs.get(db)

		db_file=None
		for file in db_info['files']:
			if file['type']=='genome':
				db_file=file['filename']

		if db_file is None:
			raise Exception(f'No genome file found for database {db}')
		
		suffixes=['amb','ann','bwt','pac','sa']
		db_files=[f"{self.db_path}/{db_file}.{x}" for x in suffixes]

		return(db_files)