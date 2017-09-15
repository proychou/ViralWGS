# HHV6 script: This script makes a new reference sequence from de novo assembled scaffolds
# Pavitra Roychoudhury
# Sep 2017

# Built to be called from hhv6_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(Biostrings);
library(RCurl);
library(parallel);

#Get latest stable version of wgs_functions.R from github
# source('./wgs_functions.R'); #or locally
script<-getURL('https://raw.githubusercontent.com/proychou/ViralWGS/master/wgs_functions.R',
							 ssl.verifypeer=FALSE)
eval(parse(text=script));

#Get args from command line 
args<-(commandArgs(TRUE));
if(length(args)==0){
	print("No arguments supplied.")
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
		print(args[[i]])
	}
}

#For testing
# bamfname<-'./contigs/ABI-HHV6A_S385_L001/aligned_scaffolds_hhv6A_ref_U1102.bam';
# reffname<-'./NC_001664.2.fasta'
# ncores<-detectCores()
	
newref<-make_ref_from_assembly(bamfname,reffname,ncores)

if(newref==FALSE) print('Failed to generate consensus from scaffolds')

