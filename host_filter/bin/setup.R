#!/bin/env Rscript

# Carried out library installation into clean conda environment which 
# just has R installed via conda

# The loadlib function will install a package if it is not already available,
# then loads it using the `library` fucnction
# If a true value is  passed via the bioc argument, the package will be
# installed from bioconductor instead of CRAN

loadlib<-function(libname, version, bioc) {

	if (!requireNamespace(libname, quietly = TRUE)) {
		if(bioc) {
		  BiocManager::install(libname)
		} else {
			install.packages(libname, version=version)
		}
 	 }
	suppressPackageStartupMessages({
		library(libname,character.only = TRUE)
	})
}
local({r <- getOption("repos")
	r["CRAN"] <- "https://cloud.r-project.org" 
	options(repos=r)
})

if (!require("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

loadlib('dplyr',   version='1.12',  bioc=FALSE)
loadlib('getopt',  version='1.20.3',bioc=FALSE)
loadlib('ggplot2', version='3.42',  bioc=FALSE)
loadlib('ggh4x', version='0.2.5',  bioc=FALSE)
loadlib('ggnewscale', version='0.4.9', bioc=FALSE)
loadlib('phyloseq',version='1.44',  bioc=TRUE)
loadlib('stringr', version='1.5.0', bioc=FALSE)
loadlib('tidyr', version='1.3.0', bioc=FALSE)
loadlib('yaml', version='2.3.7', bioc=FALSE)