#!/bin/env python

import click
import logging
from host_filter import database, workflow
from os import makedirs
from pathlib import Path

db_obj=database.database()
wf=workflow.workflow()

conf_dbs = db_obj.dbs.keys()

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
def mapping():
	wf.map()

@run.command()
@click.option('--mapq', required=True, type=click.IntRange(0,40), help='Mappinq quality threshold', 
			   default=20, show_default=True)
def filter(mapq):
	wf.filter(mapq)
	
if __name__ == "__main__":
	cli()