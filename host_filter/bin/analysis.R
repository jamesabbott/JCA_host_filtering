#!/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(operators))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(yaml))

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

optspec=matrix(c(
	'flagstat', 'f', '1', 'character', 'Path to summary table of flagstat outputs',
	'mapping',  'm', '2', 'character', 'Path to mapping file',
	'group',    'g', '2', 'character', 'Mapping file column to group samples by',
	'database', 'd', '2', 'character', 'Database name',
	'help',     'h', '0', 'logical',   'Display help'
),byrow=TRUE,ncol=5)

opt=getopt(optspec)
outdir='analysis/plots'

#TODO: Remove these!
opt$flagstat<-'analysis/mapping_summary.txt'
opt$mapping<-'mapping.txt'
opt$group<-'Microenvironment'
opt$database<-'kraken_pluspfp'

if (file.exists('config.yaml')) {
	config=read_yaml('config.yaml')
} else {
	cat('config.yaml file not found\n')
	q(status=1)
}

if (file.exists('host_filter/etc/genomes.json')) {
	genomes<-fromJSON(file='host_filter/etc/genomes.json')
	taxid<-pluck(genomes, config$database, 'taxid')
	host_species<-pluck(genomes, config$database, 'species')
} else {
	cat('etc/genomes.json file not found\n')
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
			} else {
			  metadata[[opt$group]]=factor(metadata[[opt$group]])
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

#' plot_mappings: 
#' plots proportion of mapped/unmapped read per sample
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

#' read_filtered_stats: 
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

#' plot_filtered_read_counts:
#' Generates beeswarm/boxplot per maqp across range, including unfiltered
#'
#' Available mapq data determined by presence of data within mapq_* directory
#' 
#' @param mapping_dat: dataframe of mapping data from plot_mappings

plot_filtered_read_counts<-function(mapping_dat) {
	mapq_mapping_stats<-bind_rows(map(mapqs,read_filtered_stats))

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

#' load_biom:
#' imports biom format results for a defined mapq output
#' 
#' @param db_path: path to database
#' @param mapq: mapq data set to load
#' 
#' @return biom phyloseq object

load_biom<-function(db_path, mapq) {
	biom_path<-paste0(db_path,'/biom/',mapq,'.biom')
	biom<-import_biom(biom_path)
	colnames(tax_table(biom))<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
	# clean up taxtable entries to remove k__, c__ etc prefixes
	tax_table(biom)[, colnames(tax_table(biom))] <- gsub(tax_table(biom)[, colnames(tax_table(biom))], pattern = "[a-z]__", replacement = "")

	return(biom)
}

#' load_bioms:
#' iterates through all mapq filtered outputs to read each into a ps 
#' object with `load_biom()`, resulting in a vector of ps objects. 
#' This is serialiased to an RDS file ('bioms.dat') in the database
#' directory. Should the RDS file already exist, this will be read
#' instead of the individual biom objects
#' 
#' @param database - path to database direcotry
#' 
#' @return ps-objs - vector of phyloseq objects

load_bioms<-function(database) {
	db_path=paste0((str_split_1(database,'_'))[1],'/',database)

	mapqs=list.files(db_path,pattern='mapq_[0-9]+')
	mapqs<-append(mapqs,c('unfiltered'))

	# TODO: check modtimes to reread bioms if newer than RDS file...
	if (file.exists(paste0(db_path,'/bioms.dat'))) {
		ps_objs<-readRDS(paste0(db_path,'/bioms.dat'))
	} else {
		ps_objs<-set_names(mapqs) %>% 
			map(~load_biom(.x, db_path=db_path), .progress="Loading biom files...")
		
		saveRDS(ps_objs,file=paste0(db_path,'/bioms.dat'))
	}

	return(ps_objs)
}

#' get_residual_host:
#' Determines read counts for host genome retained following filtering
#' 
#' @param mapq: mapping quality filter to process
#' @param ps_objs: vector of phyloseq objects
#' @param eukaryote: boolean to determine wether to obtain host genome 
#' or all eukaryotic sequences
#'
#' @return host_counts: dataframe of Sample, Discarded read count,
#' Retained read count, MapQ and grouping variable

get_residual_host<-function(mapq, ps_objs, eukaryote=FALSE) {

	ps<-ps_objs[[mapq]]
	
	genus<-str_split_1(host_species, ' ')[[1]]
	species<-str_split_1(host_species, ' ')[[2]]

	if (isTRUE(eukaryote) ) {
	  euks<-subset_taxa(ps, Kingdom=='Eukaryota' )
	  expr<-substitute(subset_taxa(euks, Genus != genus_val & Species != species_val),
							 list(species_val = species, genus_val=genus))
	  other_euks<-eval(expr)
	  host<-tax_glom(other_euks,taxrank = 'Kingdom')
	} else {
  	# This is evil due to how subset_taxa passes input onto subset...
	  expr <- substitute(subset_taxa(ps, Genus == genus_val & Species == species_val), 
		  list(species_val = species, genus_val = genus))
	  host<-eval(expr)
	}

	# unfiltered results need handling differently since we don't have any 
	# discarded  reads
	if (mapq=='unfiltered') {
		
		host_counts<-otu_table(host) %>%
			as.data.frame() %>%
			rename_with(str_replace, pattern='.report',replacement='') %>%
			t() %>%
			as.data.frame() %>%
			rownames_to_column('Sample') %>%
			set_colnames(c('Sample','Discarded')) %>%
			mutate('MapQ'='None') %>%
			mutate(Sample = str_replace(Sample, "_unfiltered.", "")) %>%
			left_join(metadata,join_by('Sample'=='SampleID')) 
		
	} else {

		read_counts<-read_filtered_stats(mapq) %>% 
			mutate(Total = as.numeric(Mapped)+as.numeric(Unmapped))
	
		mapped_samples<-sample_names(host)[grepl('_mapped', sample_names(host))]
		unmapped_samples<-sample_names(host)[grepl('_unmapped', sample_names(host))]
	
		expr<-substitute(subset_samples(host, Id %in% mapped_val),
					   list(mapped_val = mapped_samples))
	
		mapped_host_counts<-otu_table(eval(expr)) %>%
			as.data.frame() %>%
			rename_with(str_replace, pattern='_mapped', replacement='') %>%
			t() %>%
			as.data.frame() %>%
			rownames_to_column('Sample') %>%
			set_colnames(c('Sample','Discarded')) %>%
			mutate('MapQ'=str_replace(mapq,'mapq_','')) %>%
			mutate(Sample = str_replace(Sample, ".[0-9]*$", ""))
	
		expr<-substitute(subset_samples(host, Id %in% unmapped_val),
					list(unmapped_val = unmapped_samples))
	
		unmapped_host_counts<-otu_table(eval(expr)) %>%
			as.data.frame() %>%
			rename_with(str_replace, pattern='_unmapped',replacement='') %>%
			t() %>%
			as.data.frame() %>%
			rownames_to_column('Sample') %>%
			set_colnames(c('Sample','Retained')) %>%
			mutate('MapQ'=str_replace(mapq,'mapq_','')) %>%
			mutate(Sample = str_replace(Sample, ".[0-9]*$", ""))
	
		host_counts<-left_join(mapped_host_counts, unmapped_host_counts, join_by('Sample','MapQ')) %>%
			left_join(metadata,join_by('Sample' == 'SampleID')) 
	}
	return(host_counts)
	
}

#' plot_residual_host:
#' produces plots of host reads remaining in the samples following
#' mapq filtering, including all ps objects if provided list
#' 
#' Plots of both absolute and relative values are created
#' 
#' @param ps_objs: list of phyloseq objects
#' @param host_species: name of host organism

plot_residual_host<-function(ps_objs, host_species) {	

	mapq_labels=unlist(map(mapqs, ~str_replace(.x, 'mapq_','')))
	mapq_labels=c('None',mapq_labels)
	mapqs<-append(mapqs,'unfiltered')

	host_counts<-bind_rows(map(mapqs, ~get_residual_host(.x, ps_objs)))  %>%
		arrange(MapQ,.data[[opt$group]]) %>%
	  mutate(Type=host_species) %>%
		mutate(Sample=factor(Sample,levels=unique(Sample))) %>%
		pivot_longer(cols=c('Discarded','Retained'),names_to = 'Status', values_to = 'Count')
	
	euk_counts<-bind_rows(map(mapqs, ~get_residual_host(.x, ps_objs, eukaryote=TRUE))) %>%
		arrange(MapQ,.data[[opt$group]]) %>%
	  mutate(Type='Other Eukaryota') %>%
		mutate(Sample=factor(Sample,levels=unique(Sample))) %>%
		pivot_longer(cols=c('Discarded','Retained'),names_to = 'Status', values_to = 'Count')
	
	host_counts<-bind_rows(host_counts,euk_counts)
	
	host_counts$MapQ<-factor(host_counts$MapQ,levels=mapq_labels)
	host_counts$Type<-factor(host_counts$Type,levels=c(host_species, 'Other Eukaryota'))
	
	residual_host_plot<-ggplot(host_counts)+
		geom_beeswarm(aes(x=MapQ,y=Count,colour=Status),size=0.5, dodge.width=1, corral="wrap", corral.width=0.45)+
		facet_grid(Type ~ .data[[opt$group]]) +
		theme_bw()+
		ggtitle('Eukaryotic genome content across mapping quality range') +
		xlab('MapQ Threshold') +
		ylab('Number of Host Genome Reads') +
		theme(plot.title=element_text(size=10),
			axis.text.x=element_text(size=6),
			axis.text.y=element_text(size=6),
			axis.title=element_text(size=8),
			legend.title=element_text(size=6),
			legend.text=element_text(size=4),
			legend.position='bottom',
			strip.background = element_rect(fill='white'),
		strip.text=element_text(size=6))+ 
		scale_colour_manual(values=cbPalette)

	ggsave(residual_host_plot, file=paste0(outdir,'/residual_host.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
	
	host_counts<-host_counts %>%
		pivot_wider(names_from=Status, values_from = c('Count')) %>%
		filter(MapQ != 'None') %>%
		mutate('Proportion' = Discarded / (Discarded + Retained)) 
	
	residual_host_proportion_plot<-ggplot(host_counts)+
		geom_beeswarm(aes(x=MapQ,y=Proportion,colour=.data[[opt$group]]),size=0.5, dodge.width=1,corral="wrap", corral.width=0.5)+
		theme_bw()+
		ggtitle(paste0(str_to_title(config['database']),' genome proportion discarded across mapping quality range')) +
		xlab('MapQ Threshold') +
		ylab('Proportion') +
		theme(plot.title=element_text(size=10),
				axis.text.x=element_text(size=6),
				axis.text.y=element_text(size=6),
				axis.title=element_text(size=8),
				legend.title=element_text(size=6),
				legend.text=element_text(size=4),
				legend.position='bottom')+
		scale_colour_manual(values=cbPalette)
	
	ggsave(residual_host_proportion_plot, file=paste0(outdir,'/residual_host_proportion.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
		
}

#' Merges list of phyloseq objects into a single uber-object
#' 
#' @param ps_objs: list of phyloseq objects to merge
#' 
#' @return ps1: Merged object 

merge_ps<-function(ps_objs) {
	ps1=NULL

	for (name in names(ps_objs)) {
	  obj<-ps_objs[[name]]
	  
	  df<-sample_names(obj) %>% 
	    as.data.frame() %>%
	    rename('Sample_ID'='.') %>%
	    mutate(Filter = name) %>%
	    mutate(raw_sample_id=str_split(Sample_ID,'_',simplify=TRUE)[,1]) %>%
	    left_join(metadata,by=join_by('raw_sample_id'=='SampleID')) 
	  
	  row.names(df)=df$Sample_ID
	  
	  sample_data(obj)=df
	  
	  if (is.null(ps1)) {
	    ps1=obj
	  } else {
	    ps2=obj
	    ps1=merge_phyloseq(ps1,otu_table(ps2),sample_data(ps2),tax_table(ps2))
	  }
	}
	return(ps1)
}

#' plot_alpha
#' produces plots of alpha diversity measures. Simpson and Shannon are plotted 
#' by default. 
#' 
#' Available measures are those supported by phyloseq's plot_richness function
#' 
#' @param ps: merged phyloseq object
#' @param measures: list of measures to plot (default: Simpson and Shannon)

plot_alpha<-function(ps, measures=c('Simpson', 'Shannon')) {
	rich_plot<-plot_richness(ps, color='Filter',shape=opt$group, measures=measures)
	rich_plot<-rich_plot+
	theme_bw()+
	theme(axis.text.x=element_blank(),
		  axis.ticks.x=element_blank(),
		  plot.background = element_rect(fill='white'),
		  strip.background = element_rect(fill='white'))+
	ggtitle('Prokaryotic Alpha diversity')+
	scale_color_manual(values=cbPalette)
	ggsave(rich_plot, file=paste0(outdir,'/alpha_diversity.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
}

#' plot_beta:
#' produces plot of beta-diversity, including permanovo results
#' 
#' @param ps: merged phyloseq object

plot_beta<-function(ps) {
	ord<-ordinate(ps,"NMDS", distance='bray')

	df<-data.frame(sample_data(ps))
	bray_dist<-distance(ps, method='bray')
	res<-adonis2(bray_dist ~ Filter, data=df, permutations=5000, parallel=4)

	r2<-res$R2[[1]]
	pval<-res$`Pr(>F)`[[1]]
  
	beta_plot<-plot_ordination(ps,ord,color='Filter',shape=opt$group)
	beta_plot<-beta_plot+
	theme_bw()+
	labs(caption=paste0('PERMANOVA: r2 = ',signif(r2,3),'; p=',signif(pval,4)))

	ggsave(beta_plot, file=paste0(outdir,'/beta_diversity.pdf'), 
		device='pdf', width=15, height=10, units='cm', useDingbats=FALSE)
}

mapqs=list.files('filtered_fastq/',pattern='mapq_[0-9]+')
ps_objs<-load_bioms(opt$database)
ps<-merge_ps(ps_objs)
prokaryote_ps<-subset_taxa(ps, Kingdom == 'Bacteria' | Kingdom == 'Archaea')
mapping_dat<-plot_mappings()
plot_filtered_read_counts(mapping_dat)
plot_residual_host(ps_objs, host_species)
plot_alpha(prokaryote_ps)
plot_beta(prokaryote_ps)

