#!/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(yaml))

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")

optspec=matrix(c(
	'flagstat', 'f', '1', 'character', 'Path to summary table of flagstat outputs',
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

if (!dir.exists(outdir)) {
	dir.create(outdir)
}

plot_mappings<-function() {
	flagstat<-read.csv(opt$flagstat,sep='\t', header=TRUE) %>% 
		rename('sample'='X') 

	mapping_dat<-flagstat %>%
		mutate(mapped_prop=mapped/total) %>%
		mutate(unmapped_prop=1-mapped_prop) %>%
		select(sample,mapped_prop,unmapped_prop) %>%
		rename('Sample'='sample','Mapped'='mapped_prop','Unmapped'='unmapped_prop') %>%
		pivot_longer(cols=-c(Sample),names_to='Status')

	mapping_prop_plot<-ggplot(mapping_dat) +
		geom_bar(aes(x=Sample, y=value, fill=Status),position='fill',stat='identity')+
		ggtitle(paste('Per-sample mapping to',config['database'],'genome'))+
		ylab('Proportion')+
		theme_bw()+
 		theme(strip.background = element_rect(fill='white'),
			  legend.position = 'bottom',
			  plot.title=element_text(size=10),
			  axis.text.x=element_text(angle=45,hjust=1,size=6),
			  axis.title=element_text(size=8),
			  legend.title=element_text(size=6),
			  legend.text=element_text(size=4))+
  		scale_fill_manual(values=cbPalette)

	ggsave(mapping_prop_plot, file=paste0(outdir,'/mapping_proportion.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
}

plot_mappings()