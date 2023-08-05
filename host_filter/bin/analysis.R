#!/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(operators))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(yaml))

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

optspec=matrix(c(
	'flagstat', 'f', '1', 'character', 'Path to summary table of flagstat outputs',
	'mapping',  'm', '2', 'character', 'Path to mapping file',
	'group',    'g', '2', 'character', 'Mapping file column to group samples by',
	'help',     'h', '0', 'logical',   'Display help'
),byrow=TRUE,ncol=5)

opt=getopt(optspec)
outdir='analysis/plots'

if (file.exists('config.yaml')) {
	config=read_yaml('config.yaml')
} else {
	cat('config.yaml file not found\n')
	q(status=1)
}

if (!is.null(opt$help)) {
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$flagstat)) {
	cat('No flagstat argument provided\n')
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if ((is.null(opt$mapping) & !is.null(opt$group)) |
	((!is.null(opt$mapping) & is.null(opt$group)))) {
	cat('group (-g) and mapping (-m) arguments must both be provided to enable sample grouping\n')
	q(status=1)
} else {
	if (!is.null(opt$mapping)) {
		if (file.exists(opt$mapping)){
			metadata<-read.csv(opt$mapping,sep='\t')
			if (opt$group %!in% colnames(metadata)) {
				cat(paste0('Group ',opt$group,' not found in mapping file ',opt$mapping,'\n'))
				q(status=1)
			}
		} else {
			cat(paste0('Mapping file (',opt$mapping,') not found...\n'))
			q(status=1)
		}
	}
}

if (!dir.exists(outdir)) {
	dir.create(outdir)
}

#' Plots proportion of mapped/unmapped read per sample
#' 
#' If a mapping file and group are specified, samples will be 
#' grouped and indicated in the x-axis margin
#' 
#' @return mapping_dat

plot_mappings<-function() {
	flagstat<-read.csv(opt$flagstat,sep='\t', header=TRUE) %>% 
		rename('sample'='X') 

	mapping_dat<-flagstat %>%
		mutate(mapped_prop=mapped/total) %>%
		mutate(unmapped_prop=1-mapped_prop) %>%
		select(sample, mapped_prop, unmapped_prop) %>%
		rename('Sample'='sample', 'Mapped'='mapped_prop','Unmapped'='unmapped_prop') 

	if (!is.null(opt$group)) {
	 	mapping_dat <- mapping_dat %>% 
	 		left_join(metadata,join_by('Sample'=='SampleID')) %>%
			pivot_longer(cols=-c('Sample', opt$group),names_to='Status') %>%
			arrange(.data[[opt$group]]) %>%
			mutate(Sample=factor(Sample,levels=unique(Sample)))
	} else {
		mapping_dat<-mapping_dat %>%
			pivot_longer(cols=-c('Sample'),names_to='Status')
	}
	
	mapping_prop_plot<-ggplot(mapping_dat) +
		geom_bar(aes(x=Sample, y=value, fill=Status),position='fill',stat='identity')+
		ggtitle(paste('Per-sample mapping to',config['database'],'genome'))+
		ylab('Proportion')+
		theme_bw()+
 		theme(legend.position = 'bottom',
			  plot.title=element_text(size=10),
			  axis.text.x=element_text(angle=45,hjust=1,size=6),
			  axis.text.y=element_text(size=6),
			  axis.title=element_text(size=8),
			  legend.title=element_text(size=6),
			  legend.text=element_text(size=4))+
		scale_y_continuous(limits = c(0, 1)) +
		scale_fill_manual(values=cbPalette[1:2])

	if (!is.null(opt$group)) {
		mapping_prop_plot<-mapping_prop_plot+
		new_scale('fill')+
		scale_fill_manual(values=cbPalette[3:8])+
		geom_tilemargin(aes(x=Sample, y=0, fill=.data[[opt$group]]),outside=FALSE,sides='b')+
		coord_cartesian(clip = "off")	
	}
	ggsave(mapping_prop_plot, file=paste0(outdir,'/mapping_proportion.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
	
	return(mapping_dat)
}

#' Reads per-sample mapping counts for a given mapq directory
#' 
#' @param mapq: mapq directory name
#' 
#' @return DataFrame, with columns Sample, Mapped, Unmapped and Mapq

read_filtered_stats<-function(mapq) {
	files<-list.files(paste0('filtered_fastq/',mapq), pattern=".*_stats.txt")
	stats<-map(files, function(file) {
		sample_id<-str_replace(file,'_stats.txt','')
		dat<-as.data.frame(t(read.table(paste0('filtered_fastq/',mapq,'/',file), sep="\t",header=FALSE)))

		colnames(dat)<-dat[1,]
		dat<-dat[-1,]
		dat['Mapped']=str_trim(dat['Mapped'])
		dat['Sample']<-sample_id
		dat
	})

	results<-bind_rows(stats)
	results$MapQ<-mapq

	return(results)
}

plot_filtered_read_counts<-function(mapping_dat) {
	mapqs=list.files('filtered_fastq/',pattern='mapq_[0-9]+')
	mapq_mapping_stats<-bind_rows(map(mapqs,read_filtered_stats))# %>%

	mapq_mapping_stats<-mapq_mapping_stats %>%
		mutate(mapped_prop=as.numeric(Mapped)/(as.numeric(Mapped)+as.numeric(Unmapped))) %>%
		mutate(MapQ=str_replace(MapQ,'mapq_','')) %>%
		left_join(metadata,join_by('Sample'=='SampleID')) %>%
		select(c(Sample,MapQ,mapped_prop,opt$group))
	
	# Add data for non-mapq filtered reads...
	unfiltered_stats<-mapping_dat %>%
		filter(Status=='Mapped') %>%
		mutate(MapQ='None') %>%
		mutate(value=as.numeric(value)) %>%
		rename('mapped_prop'='value') 

	# Setup a factor to order mapq with 'None' first...
	mapq_mapping_stats<-bind_rows(mapq_mapping_stats,unfiltered_stats)
	mapq_vals<-unique(mapq_mapping_stats$MapQ) 
	mapq_vals<-c('None',mapq_vals[mapq_vals %!in% c('None')])
	mapq_mapping_stats$MapQ<-factor(mapq_mapping_stats$MapQ, levels=mapq_vals)

	mapq_mapping_plot<-ggplot(mapq_mapping_stats)+
		geom_beeswarm(aes(x=MapQ,y=mapped_prop,colour=.data[[opt$group]]),dodge.width=0.7) +
		geom_boxplot(aes(x=MapQ,y=mapped_prop,colour=.data[[opt$group]]),position=position_dodge2(width=0.7),alpha=0.1,linewidth=0.2,notch=TRUE) +
		theme_bw() +
		ggtitle('Proportion of discarded reads across mapping quality range') +
		xlab('MapQ Threshold') +
		ylab('Discarded Proportion of Reads') +
		theme(plot.title=element_text(size=10),
			  axis.text.x=element_text(size=6),
			  axis.text.y=element_text(size=6),
			  axis.title=element_text(size=8),
			  legend.title=element_text(size=6),
			  legend.text=element_text(size=4))+
		scale_colour_manual(values=cbPalette)

	ggsave(mapq_mapping_plot, file=paste0(outdir,'/mapq_mapping_proportion.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
}

mapping_dat<-plot_mappings()
plot_filtered_read_counts(mapping_dat)