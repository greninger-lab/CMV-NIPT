library(svMisc)
library(seqinr)
library(Rsamtools)
library(ggplot2)

setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download')

blast_positive<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/blast_positive_sequence_info.csv')
all_samples<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/all_sample_data.csv')


intermediates<-all_samples[all_samples$classification=='intermediate positive',]

blast_positive$sample<-NA
for(i in 1:nrow(blast_positive)){ 
  progress(i)  
  blast_positive$sample[i]<-strsplit(as.character(blast_positive$nameslist[i]),'[.]')[[1]][1]
}

intermediates_with_sequences<-setNames(data.frame(matrix(ncol = ncol(blast_positive), nrow = 0)), colnames(blast_positive))
for(i in 1:nrow(intermediates)){ 
  temp<-which(intermediates$sample[i]==blast_positive$sample)
  print(temp)
  for(j in 1:length(temp)){ 
    intermediates_with_sequences<-rbind(intermediates_with_sequences, blast_positive[temp[j],])
    
    }
}
intermediates_with_sequences$nameslist<-as.character(intermediates_with_sequences$nameslist)
for(i in 1:nrow(intermediates_with_sequences)){ 
  #print(as.character(strsplit(as.character(intermediates_with_sequences$nameslist[i]),'-')[[1]][2]))
  intermediates_with_sequences$nameslist[i]<-as.character(strsplit(as.character(intermediates_with_sequences$nameslist[i]),'-')[[1]][2])
  }

write.csv(intermediates_with_sequences,'intermediate_positive_sequences.csv')

system('mkdir intermediate_positive_fastas')
setwd('intermediate_positive_fastas/')

for(i in 1:nrow(intermediates_with_sequences)){ 
  temp_header<-paste0(intermediates_with_sequences$sample[i],intermediates_with_sequences$nameslist[i])
  temp_seq<-intermediates_with_sequences$seqslist[i]
  print(temp_header)
  temp_file_name<-paste(i,'.tempfa')
  write.fasta(sequences = temp_seq, names = temp_header, file.out = temp_file_name)
}
system('cat *.tempfa > intermediate_sequences.fasta')
system('rm *.tempfa')

intermediates_with_sequences$isize<-NA
for(i in 1:nrow(intermediates_with_sequences)){ 
  #progress(i,nrow(intermediates_with_sequences))
  tempbamname<-paste0('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/original_bams/',intermediates_with_sequences$sample[i],'.sam.bam')
  tempbam<-scanBam(tempbamname)
  read<-which(grepl(intermediates_with_sequences$nameslist[i], tempbam[[1]]$qname))
  temp_isize<-((tempbam[[1]]$isize[read]))
  if( any(is.na(temp_isize))) {
    temp_isize<-0
  }
  temp_isize<-abs(temp_isize[1])
  print(temp_isize)
  intermediates_with_sequences$isize[i]<-temp_isize
  }
intermediates_median<-median(intermediates_with_sequences$isize)
ggplot(intermediates_with_sequences, aes(x = intermediates_with_sequences$isize)) + 
  geom_freqpoly() + 
  theme_classic() 
