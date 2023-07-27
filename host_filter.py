#!/bin/env python

import click
import logging
from glob import glob
from host_filter import genome, workflow, classify, analysis
from os import makedirs
from pathlib import Path
import re

db_obj=genome.database()
wf=workflow.workflow()
classify_obj=classify.classify()
analysis_obj=analysis.analysis() 

conf_dbs = db_obj.dbs.keys()
conf_ref_dbs = classify_obj.dbs.keys()

def validate_range(ctx,param,value):
	match=re.search(r'^([0-9]+)-([0-9]+):([0-9]+)$',value)
	if match:
		return(match.group(1),match.group(2),match.group(3))
	else:
		raise click.BadParameter("Range format must be start-end:increment i.e. 10-30:5")
		
		
@click.group()
def cli():
	makedirs('logs',exist_ok=True)
	logname=Path(f'logs/host_filter.log')

	logging.basicConfig(
		filename=logname,
		filemode='a',
		format='%(asctime)s, %(name)s %(levelname)s %(message)s',
		datefmt='%H:%M:%S',
		encoding='UTF-8',
		level=logging.DEBUG
	)
	global logger
	logger = logging.getLogger(__name__)

@cli.group()
def genome():
	pass

@genome.command()
def status():
	db_obj.available_dbs()

@genome.command()
@click.option('--db', required=True, type=click.Choice(conf_dbs), help='database name')
def download(db):
	db_obj.download(db)
	db_obj.bwa_index(db)

@genome.command()
@click.option('--db', required=True, type=click.Choice(conf_dbs), help='database name')
def index(db):
	db_obj.bwa_index(db)

@genome.command()
@click.option('--db', required=True, type=click.Choice(conf_dbs), help='database name')
def clean(db):
	db_obj.clean(db)

@cli.group()
def db():
	pass

@db.command()
def status():
	classify_obj.available_dbs()

@db.command()
@click.option('--db', required=True, type=click.Choice(conf_ref_dbs), help='reference database name')
def download(db):
	classify_obj.download(db)

@cli.group()
def run():
	pass

@run.command()
@click.option('-v','--verbose', is_flag=True,help='Generate verbose output')
def mapping(verbose):
	wf.map(verbose)

@run.command()
@click.option('-m','--mapq', required=True, type=click.IntRange(0,40), help='Mappinq quality threshold', 
			   default=20, show_default=True)
@click.option('-v','--verbose', is_flag=True,help='Generate verbose output')
def filter(mapq,verbose):
	wf.filter(mapq,verbose)

@run.command()
@click.option('-m','--mapq_range', required=True, callback=validate_range,
	      help='Range of mapping qualities to filter i.e. "10-30:5" will carry out separate filtering steps between MAPQ 10 and 30 with a increment of 5')
@click.option('-v','--verbose', is_flag=True,help='Generate verbose output')
def range_filter(mapq_range,verbose):
	start_mapq, end_mapq, increment=mapq_range
	#increment end_mapq by 1 to ensure the final value is included if it on a boundary
	wf.range_filter(int(start_mapq),int(end_mapq)+1, int(increment), verbose)

@run.command()
@click.option('-v','--verbose', is_flag=True,help='Generate verbose output')
def metaphlan(verbose):
	wf.metaphlan(verbose)

@run.command()
@click.option('-v','--verbose', is_flag=True, help='Generate verbose output')
@click.option('-d','--db', type=click.Choice(conf_ref_dbs), default='kraken_standard', help='reference database name')
def kraken(db,verbose):
	wf.kraken(db,verbose)

@cli.command()
@click.option('-r','--results', type=click.Choice(conf_ref_dbs), default='kraken_standard', help='result set to analyse')
@click.option('-v','--verbose', is_flag=True, help='Generate verbose output')
def analyse(results,  verbose):
	results_path=f"{results.split('_')[0]}/{results}"
	mapq_range=glob(f"{results_path}/*")
	mapq_range=[Path(x).name for x in mapq_range]
	for mapq in mapq_range:
		analysis_obj.analyse(results,mapq,verbose)

if __name__ == "__main__":
	cli()