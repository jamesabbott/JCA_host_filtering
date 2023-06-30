#!/bin/env python

import click
import logging
from host_filter import database, workflow
from os import makedirs
from pathlib import Path
import re

db_obj=database.database()
wf=workflow.workflow()

conf_dbs = db_obj.dbs.keys()

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
def db():
	pass

@db.command()
def status():
	db_obj.available_dbs()

@db.command()
@click.option('--db', required=True, type=click.Choice(conf_dbs), help='database name')
def download(db):
	db_obj.download(db)
	db_obj.bwa_index(db)

@db.command()
@click.option('--db', required=True, type=click.Choice(conf_dbs), help='database name')
def index(db):
	db_obj.bwa_index(db)

@db.command()
@click.option('--db', required=True, type=click.Choice(conf_dbs), help='database name')
def clean(db):
	db_obj.clean(db)

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
@click.option('-r','--range', required=True, callback=validate_range,
	      help='Range of mapping qualities to filter i.e. "10-30:5" will carry out separate filtering steps between MAPQ 10 and 30 with a increment of 5')
@click.option('-v','--verbose', is_flag=True,help='Generate verbose output')
def range_filter(range,verbose):
	start_mapq, end_mapq, increment=range
	wf.range_filter(int(start_mapq),int(end_mapq), int(increment), verbose)

if __name__ == "__main__":
	cli()