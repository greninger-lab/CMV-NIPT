library(svMisc)

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
