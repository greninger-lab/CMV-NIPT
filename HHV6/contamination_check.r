library(svMisc)
library(seqinr)
library(Rsamtools)
library(ggplot2)

setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download')

blast_positive<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/blast_positive_sequence_info.csv')
all_samples<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/all_sample_data.csv')

blast_positive$sample<-NA
for(i in 1:nrow(blast_positive)){ 
  progress(i)  
  blast_positive$sample[i]<-strsplit(as.character(blast_positive$nameslist[i]),'[.]')[[1]][1]
}


blast_positive$nameslist<-as.character(blast_positive$nameslist)
for(i in 1:nrow(blast_positive)){ 
  #print(as.character(strsplit(as.character(intermediates_with_sequences$nameslist[i]),'-')[[1]][2]))
  blast_positive$nameslist[i]<-as.character(strsplit(as.character(blast_positive$nameslist[i]),'-')[[1]][2])
}

#filling in insert sizes
blast_positive$pos<-NA
for(i in 1:nrow(blast_positive)){ 
  #progress(i,nrow(intermediates_with_sequences))
  tempbamname<-paste0('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/original_bams_deduplicated/',blast_positive$sample[i],'.sam.bam.sorted.bam_dedup.bam')
  tempbam<-scanBam(tempbamname)
  read<-which(grepl(blast_positive$nameslist[i], tempbam[[1]]$qname))
  temp_pos<-((tempbam[[1]]$pos[read]))
  if( any(is.na(temp_pos))) {
    temp_pos<-0
  }
  temp_pos<-abs(temp_pos[1])
  print(temp_pos)
  blast_positive$pos[i]<-temp_pos
}
blast_positive<-blast_positive[complete.cases(blast_positive$pos),]


