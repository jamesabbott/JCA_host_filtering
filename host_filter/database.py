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

def read_dbs():
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

def available_dbs():
	"""
	Lists databases available for download

	Required parameters: 
		None

	Returns: 
		None
	"""
	dbs=read_dbs()
	output=[['Database','Downloaded','BWA Indexed']]
	for db in dbs:
		db_info=dbs[db]
		count=0
		found=0
		bwa_index='No'
		have_db='No'

		for file in db_info['files']:
			count+=1
			if (Path('databases')/Path(file['filename']+'.download_complete')).exists():
				found+=1
			if (Path('databases')/Path(file['filename']+'.sa')).exists():
				bwa_index='Yes'

		if found==count:
			have_db='Yes'
		
		output.append([db,have_db,bwa_index])
		
	print(tabulate(output,headers='firstrow',tablefmt='fancy_grid'))

def download(db):
	"""
	Downloads files for specified database 

	Required parameterse:
		db(str): database name

	Returns: 
		None
	"""
	dbs=read_dbs()
	db_info=dbs.get(db)
	makedirs(Path('databases'),exist_ok=True)
	with yaspin(text=f"Downloading {db} database...", timer=True):
		for file in db_info['files']:

			dl_file=Path('databases') / Path(file['filename']+'.dl')
			final_file=Path('databases')/Path(file['filename'])
			done_file=Path('databases') / Path(file['filename']+'.download_complete')

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

def bwa_index(db):
	"""
	Generates a BWA index for the specified database

	Required params: 
		db(str): database name
	
	Returns:
		None
	"""

	dbs=read_dbs()
	db_info=dbs.get(db)
	with yaspin(text=f"BWA indexing {db} database...", timer=True):
		for file in db_info['files']:
			if file['type']=='genome':
				filename=file['filename']
				break

		genome_file=Path('databases')/Path(filename)
		cmd=f'bwa index {genome_file}'
		common.run_command(cmd,f'bwa_index_{db}')

def clean(db):
	"""
	Removes files relating to specified database

	Required parameters:
		db(str): database name

	Returns:
		None
	"""
	dbs=read_dbs()
	db_info=dbs.get(db)
	with yaspin(text=f"Cleaning up {db} database..."):
		logging.info(f'{db}: Database removed...')
		for file in db_info['files']:
			search_path=Path('databases')/Path(file['filename']+'*')
			db_files=glob(str(search_path))
			for db_file in db_files:
				remove(db_file)