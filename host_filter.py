#!/bin/env python

import click
from host_filter import database

@click.group()
def cli():
	pass

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
	
if __name__ == "__main__":
	cli()