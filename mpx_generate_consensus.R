# HSV script: This script imports bam files and makes a consensus sequence
# Pavitra Roychoudhury
# Oct 2017

# Built to be called from hhv6_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(ShortRead);
library(Biostrings);
library(RCurl);

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

#For testing (these args should come from command line)
# s1<-'/fh/fast/jerome_k/HSV_WGS/fastq_files/NGS09_HSV1_morereads/2016-01040_S451_L001_R1_001.fastq'

#Files, directories, target site
merged_bam_folder<-'./remapped_reads/'; 
mapped_reads_folder<-'./mapped_reads/';
con_seqs_dir<-'./consensus_seqs_all';

#Make consensus sequences against reference--returns TRUE if this worked
conseq<-clean_consensus_mpx(sampname,merged_bam_folder,mapped_reads_folder,ref);

#Prepare seqs for annotation 
if(conseq==TRUE){
  if(!dir.exists('./annotations_prokka')) dir.create('./annotations_prokka');
  
  #Remove all Ns at the beginning and end of the seq, write to folder
  fname<-grep(sampname,list.files(con_seqs_dir,full.names=T),value=T);
  con_seq<-readDNAStringSet(fname);
  con_seq_trimmed<-DNAStringSet(gsub("N*N$",'',gsub("^N*",'',as.character(con_seq))));
  names(con_seq_trimmed)<-substring(names(con_seq),1,20); #prokka needs contig name to be <=20 chars long
  sampdir<-paste('./annotations_prokka/',sampname,sep='');
  if(!dir.exists(sampdir)) dir.create(sampdir); #create folder for the sample
  writeXStringSet(con_seq_trimmed,file=paste(sampdir,'/',sampname,'.fa',sep=''),format='fasta');
  
}else{
  print('Failed to generate consensus sequences.')
}

