#!/bin/env python

import click
import logging
from host_filter import database,workflow
from os import makedirs
from pathlib import Path
import snakemake

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
	database.available_dbs()

@db.command()
@click.option('--db', help='database name')
def download(db):
	database.download(db)
	database.bwa_index(db)

@db.command()
@click.option('--db', help='database name')
def index(db):
	database.bwa_index(db)

@db.command()
@click.option('--db', help='database name')
def clean(db):
	database.clean(db)

@cli.command()
@click.argument('action', type=click.Choice(['mappings','filter']))
def run(action):
	workflow.run(action)
	
if __name__ == "__main__":
	cli()