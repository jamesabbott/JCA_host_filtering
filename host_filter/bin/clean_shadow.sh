#!/bin/env bash

hosts=$(qhost|tail -n+4 | awk '{print $1}')

for host in ${hosts[@]}; do
	status=$(qhost -h ${host}|tail -n+4|awk '{print $7}')
	echo $status
	if [[ "$status" != "-" ]]; then
		qsub -j y -mods l_hard hostname $host -b y rm -rfv /tmp/shadow
	else
		echo "$host is down"
	fi

done