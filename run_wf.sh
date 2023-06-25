#!/bin/env bash
snakemake --snakefile host_filter/etc/Snakefile -r --rerun-incomplete --use-conda \
	--cluster-config cluster_config.json --configfile config.yaml -k \
	--shadow-prefix=/tmp --latency-wait 60 -j 100 --jn host_filter_{jobid} --drmaa \
	" -cwd -j y -o logs -jc {cluster.jc} -pe smp {cluster.threads} -mods l_hard mfree {cluster.mem_free} -adds l_hard local_free {cluster.local_free} -adds l_hard ramdisk {cluster.ramdisk} " 

