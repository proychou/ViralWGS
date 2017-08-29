#Collection of functions for working with WGS data
#Pavitra Roychoudhury
#Aug 2017

#Return the number of mapped reads in a bam file
n_mapped_reads<-function(bamfname){
  require(Rsamtools)
  indexBam(bamfname)
  if(file.exists(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    return(idxstatsBam(bamfname)$mapped)
  }else{
    return(NA)
  }
}

#Takes in a bam file, produces consensus sequence
generate_consensus<-function(bamfname){
  if(!is.na(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    
  	#Index bam if required
    if(!file.exists(paste(bamfname,'.bai',sep=''))){
      baifname<-indexBam(bamfname); 
    }else{
      baifname<-paste(bamfname,'.bai',sep='');
    }
    
    #Import bam file
    params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                         what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
    gal<-readGAlignments(bamfname,index=baifname,param=params);
    # summary(gal);
    
    #Remove any contigs with mapq <2
    # gal<-gal[mcols(gal)$mapq>2];
    
    #Fill gaps with Ns and generate consensus
    qseq_on_ref<-sequenceLayer(mcols(gal)$seq,cigar(gal),from="query",to="reference");
    cm<-consensusMatrix(qseq_on_ref,as.prob=T,shift=start(gal)-1,width=seqlengths(gal))
    cm["+", colSums(cm) == 0] <- 1
    con_seq<-DNAStringSet(consensusString(cm, ambiguityMap='N'));
    con_seq<-DNAStringSet(gsub('\\+','N',con_seq));
    names(con_seq)<-sub('.bam','_consensus',basename(bamfname));
      
    #Delete bai file
    file.remove(baifname);
    
    return(con_seq);
  }else{
    return(NA)
  }
}

clean_consensus_hsv<-function(sampname,merged_bam_folder){
  require(Rsamtools); 
  require(GenomicAlignments);
  require(Biostrings);
  mapping_stats<-data.frame(ref=c('hsv1_ref','hsv2_sd90e','hsv2_ref_hg52'),
                            bamfname_merged=c(grep(sampname,list.files(merged_bam_folder,'_hsv1_ref*.bam$',full.names=T),value=T),
                                              grep(sampname,list.files(merged_bam_folder,'_hsv2_sd90e*.bam$',full.names=T),value=T),
                                              grep(sampname,list.files(merged_bam_folder,'_hsv2_ref_hg52*.bam$',full.names=T),value=T)),
                            perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);
  
  #Import mapped reads + assembly and generate consensus
  con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
  if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
  dummyvar<-lapply(con_seqs,function(x)
    writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
  rm(dummyvar)
  
  #Compute %Ns
  mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
  mapping_stats$width<-unlist(lapply(con_seqs,width));
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  if(!dir.exists('./stats/')) dir.create('./stats/');
  write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
  
  #Write the best consensus
  if(!dir.exists('./cleaned_consensus')) dir.create('./cleaned_consensus');
  con_seq<-con_seqs[[which.min(mapping_stats$perc_Ns)]];
  writeXStringSet(con_seq,file=paste('./cleaned_consensus/',names(con_seq),'.fasta',sep=''),format='fasta');
  
  return(con_seq)
}

clean_consensus_hhv6<-function(sampname,merged_bam_folder,mapped_reads_folder){
	require(Rsamtools); 
	require(GenomicAlignments);
	require(Biostrings);
	mapping_stats<-data.frame(ref=c('hhv6A_ref_U1102','hhv6B_ref_z29'),
														bamfname_merged=c(grep(sampname,list.files(merged_bam_folder,'_hhv6A_ref_U1102*.bam$',full.names=T),value=T),
																							grep(sampname,list.files(merged_bam_folder,'_hhv6B_ref_z29*.bam$',full.names=T),value=T)),
														bamfname_mapped=c(grep(sampname,list.files(mapped_reads_folder,'_hhv6A_ref_U1102.*bam$',full.names=T),value=T),
																							grep(sampname,list.files(mapped_reads_folder,'_hhv6B_ref_z29.*bam$',full.names=T),value=T)),
														mapped_reads=0,perc_Ns=0,num_Ns=0,width=0,
														stringsAsFactors=F);
	
	#Import mapped reads + assembly and generate consensus
	con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
	if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
	dummyvar<-lapply(con_seqs,function(x)
		writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
	rm(dummyvar)
	
	#Compute #mapped reads and %Ns
	mapping_stats$mapped_reads<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
	mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
	mapping_stats$width<-unlist(lapply(con_seqs,width));
	mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
	if(!dir.exists('./stats/')) dir.create('./stats/');
	write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
	
	return(TRUE)
}


#Compute coverage stats
get_coverage<-function(bamfname){
  require(Rsamtools);
  require(GenomicAlignments);
  
  if(file.exists(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    #Import alignment
    if(file.exists(paste(bamfname,'.bai',sep='')))
      file.remove(paste(bamfname,'.bai',sep='')); #remove any old index files
    baifname<-indexBam(bamfname); #Make an index file
    params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                         what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
    gal<-readGAlignments(bamfname,index=baifname,param=params);
    # summary(gal);
    cov<-coverage(gal);
    mapped<-length(gal);
    avg_cov<-mean(cov);
    sd_cov<-sd(cov);
    min_cov<-min(cov);
    max_cov<-max(cov);
    file.remove(baifname);
  }else{
    mapped<-NA;
    avg_cov<-NA;
    sd_cov<-NA;
    min_cov<-NA;
    max_cov<-NA;
  }
  return(data.frame(mapped,avg_cov,sd_cov,min_cov,max_cov))
}

#Extracts number of reads and read widths from html report generated by fastqc
fastqc_readstats<-function(fname){
	require(rvest)
  if(file.exists(fname)){
    tmp_fastqc<-read_html(fname);
    tmp_table<-html_table(tmp_fastqc)[[1]];
    fastq_reads<-as.numeric(tmp_table[tmp_table$Measure=='Total Sequences','Value']);
    fastq_width<-tmp_table[tmp_table$Measure=='Sequence length','Value']; #returns single number for raw reads and range for trimmed
    gc<-as.numeric(tmp_table[tmp_table$Measure=='%GC','Value']);
  }else{
    fastq_reads<-NA;
    fastq_width<-NA;
    gc<-NA;
  }
  return(data.frame(fastq_reads,fastq_width,gc,stringsAsFactors=F));
}


#Compute stats on a consensus seq (or really any fasta file)
conseq_stats<-function(fname){
  if(file.exists(fname)){
    conseq<-readDNAStringSet(fname,format='fasta');
    width<-width(conseq);
    Ns<-sum(letterFrequency(conseq,c('N','+')));
    percNs<-100*Ns/width;
  }else{
    width<-NA; Ns<-NA; percNs<-NA;
  }
  return(data.frame(width,Ns,percNs));
}



#VCF to data frame for a vcf generated by Lofreq
vcf_to_df<-function(vcf_fname,sampid){
	require(VariantAnnotation);
	vcf<-readVcf(vcf_fname);
	results<-data.frame(samp_id=sampid,pos=start(rowRanges(vcf)),af=info(vcf)$AF,dp=info(vcf)$DP,ref=ref(vcf),
											alt=unlist(alt(vcf)),stringsAsFactors=F);
	results$snpid<-paste(results$ref,'_',results$pos,'_',results$alt,sep='');
	results$major_af<-unlist(lapply(results$af,function(x)max(x,1-x)));
	results$minor_af<-unlist(lapply(results$af,function(x)min(x,1-x)));
	return(results)
}

