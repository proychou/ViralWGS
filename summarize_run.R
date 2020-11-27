#Run stats
library(parallel); library(RCurl); library(lubridate); library(rvest)
library(Rsamtools);library(GenomicAlignments); library(tidyverse)

#Get latest stable version of wgs_functions.R from github
# source('./wgs_functions.R'); #or locally
script<-getURL('https://raw.githubusercontent.com/proychou/ViralWGS/master/wgs_functions.R',
               ssl.verifypeer=FALSE)
eval(parse(text=script));

#args from command line should include rundir, fqdir, ncores
args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
    print(args[[i]])
  }
}

#list all fastq files in the run
all_fastq_files<-list.files(fqdir,pattern='.fastq.gz',full.names=T,recursive=F);
length(all_fastq_files);

#only R1s to match names
all_samples<-data.frame(fname_full=grep("_R1_|_R1.fastq.gz|_1.fastq.gz",
                                        all_fastq_files,value=T),stringsAsFactors=F);
nrow(all_samples);

#extract filenames and sample ids
all_samples$fname_short<-unlist(lapply(basename(all_samples$fname_full),function(x)
  strsplit(x,'_R1_*|_1.fastq*')[[1]][1]));
head(all_samples);

#Columns of database
cols<-c('SpID','seq_run','pairedend','fastq_reads','fastq_width',
        'mapped_reads','avg_cov','max_cov','conNsperc','rundate')
all_samples[,cols]<-0;
all_samples$SpID<-str_split(all_samples$fname_short,'_',simplify=T)[,1]
all_samples$seq_run<-str_split(
  str_split(all_samples$fname_full,'fastq_files/',simplify=T)[,2],'/',simplify=T)[,1]


for(i in 1:nrow(all_samples)){
  print(paste('File',i,'of',nrow(all_samples),':',all_samples$SpID[i],'...'));  
  
  #Check paired
  if(file.exists(sub('_R1','_R2',all_samples$fname_full[i]))){
    all_samples$pairedend[i]<-TRUE 
  }else all_samples$pairedend[i]<-FALSE;
  
  
  #Fastq info from fastqc stats
  fqc_fnames<-grep(paste0(all_samples$fname_short[i],'_'),
                   list.files(paste0(rundir,'/fastqc_reports_raw'),
                              pattern='*_fastqc.html',full.names=T),value=T);
  tmp_table<-do.call('rbind',mclapply(fqc_fnames,function(x)
    fastqc_readstats(x),mc.cores=ncores));
  if(!is.null(tmp_table)){
    all_samples$fastq_reads[i]<-sum(tmp_table$fastq_reads);
    all_samples$fastq_width[i]<-paste(tmp_table$fastq_width,collapse=';')
  }else{
    all_samples$fastq_reads[i]<-NA;
    all_samples$fastq_width[i]<-NA;
  }
  
  #Coverage 
  bamfname<-grep(paste0(strsplit(all_samples$fname_short[i],'_')[[1]][1],'.bam'),
                 list.files(paste0(rundir,'/mapped_reads'),
                            '*.bam$',full.names=T),value=T);
  if(length(bamfname)>0){
    cov<-get_coverage(bamfname);
    all_samples[i,c('mapped_reads','avg_cov','max_cov')]<-
      cov[,c('mapped','avg_cov','max_cov')];
    rm(cov)
  }else if(length(bamfname)==0){
    all_samples[i,c('mapped_reads','avg_cov','max_cov')]<-NA;
  }else{
    all_samples[i,c('mapped_reads','avg_cov','max_cov')]<-'duplicate SpID in run--check!';
  }; rm(bamfname)
  
  
  #Ns in consensus
  csfname<-grep(paste0(strsplit(all_samples$fname_short[i],'_')[[1]][1],'.fasta'),
                list.files(paste0(rundir,'/consensus_seqs'),'*.fasta',
                           full.names=T),value=T);
  if(length(csfname)>0){
    cstats<-conseq_stats(csfname);
    all_samples$conNsperc[i]<-cstats$percNs;
  }else if(length(csfname)==0){
    all_samples$conNsperc[i]<-NA;
  }else{
    all_samples$conNsperc[i]<-'duplicate SpID in run--check!';
  }
  
  
  #analysis run date: use timestamp on consensus seq
  if(any(is.na(all_samples[i,]))){ #Determine whether the run failed or if there were too few reads
    if(!is.na(all_samples$fastq_reads[i])){
      if(all_samples$fastq_reads[i]<5000){
        all_samples$rundate[i]<-'Too few reads'
      }else if(is.na(all_samples$mapped_reads[i])){
        all_samples$rundate[i]<-'Run failed'
      }else if(all_samples$avg_cov[i]<10){
        all_samples$rundate[i]<-'Low coverage'
      }else if(is.na(all_samples$conNsperc[i])){
        all_samples$rundate[i]<-'No consensus seq'
      }else if(all_samples$mapped_reads[i]<5000){
        all_samples$rundate[i]<-'Low on tgt reads'
      }
    }else{
      all_samples$rundate[i]<-'Run failed'
    }
    
  }else{
    all_samples$rundate[i]<-as.character(file.info(grep(
      all_samples$SpID[i],list.files(
        paste0(rundir,'/consensus_seqs'),'*.fasta',full.names=T),
      value=T))$mtime)
  }
  
  #Write to file just in case it crashes before done
  write_csv(all_samples,'./RunSummary.csv');
    
}; rm(i);
write_csv(all_samples,'./RunSummary.csv');



