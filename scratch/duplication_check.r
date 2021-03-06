library(ggplot2)
library(Rsamtools)
library(svMisc)


ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

resequenced_verified_reads<-ReadFasta('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/resequenced_repeatmasked_7_removed.fasta')

#coverage graph: 

resequenced_bam<-scanBam('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/resequenced.bam.filtered.sam.bam')

resequenced_verified_reads$name_trimmed<-NA
for(i in 1:nrow(resequenced_verified_reads)){ 
  progress(i,nrow(resequenced_verified_reads))
  resequenced_verified_reads$name_trimmed[i]<-strsplit(as.character(resequenced_verified_reads$name[i]),'/')[[1]][1]
}

#shouldn't be making this unique values, pairs count as coverage 
to_remove<- which(duplicated(as.character(resequenced_verified_reads$name_trimmed)))
unique_resequenced_verified_reads<-resequenced_verified_reads[-to_remove,]

#resequenced_bam_names<-resequenced_bam[[1]]$qname
# resequenced_bam_names_trimmed<-c()
# for(i in 1:length(resequenced_bam_names)){ 
#   progress(i,nrow(verified_reads))
#   resequenced_bam_names_trimmed <-append(resequenced_bam_names_trimmed,strsplit(as.character(resequenced_bam_names[i]),'/')[[1]][1])
# }
# 

resequenced_verified_reads$pos<-NA


unique_resequenced_verified_reads$paired<-NA
for(i in 1:nrow(unique_resequenced_verified_reads)){ 
  tempindex<-which(as.character(unique_resequenced_verified_reads$name_trimmed[i]) == resequenced_bam[[1]]$qname)
  if(is.na(tempindex[2])){ 
    #print(index)
    unique_resequenced_verified_reads$paired[i]<-FALSE
  }
  else{
    print(tempindex)
    unique_resequenced_verified_reads$paired[i]<-TRUE
  }
}





positions_list<-c()
read_id_list<-c()
for(i in 1:length(unique_resequenced_verified_reads$name_trimmed)){
  index<- which(as.character(unique_resequenced_verified_reads$name_trimmed[i]) == resequenced_bam[[1]]$qname)
  if(unique_resequenced_verified_reads$paired[i]==FALSE){ 
    #  print((index))
    positions_list<-append(positions_list,resequenced_bam[[1]]$pos[index])
    read_id_list<-append(read_id_list,resequenced_bam[[1]]$qname[index])
  }
  else{
    #print(index)
    # if(resequenced_bam[[1]]$isize[index[1]] > 0 ){ 
    positions_list<-append(positions_list,resequenced_bam[[1]]$pos[index[1]])
    read_id_list<-append(read_id_list,resequenced_bam[[1]]$qname[index[1]])
    
    positions_list<-append(positions_list,resequenced_bam[[1]]$pos[index[2]])
    read_id_list<-append(read_id_list,resequenced_bam[[1]]$qname[index[2]])
    
    #}
  }
  
}

positions_with_read_ids<-data.frame(positions_list,read_id_list)


positions_with_read_ids$paired<-NA
for(i in 1:nrow(positions_with_read_ids)){ 
  #progress(i,nrow(positions_with_read_ids))
  progress(i,nrow(positions_with_read_ids))
  tempindex<-which(grepl(as.character(positions_with_read_ids$read_id_list[i]), as.character(unique_resequenced_verified_reads$name_trimmed)))
  #print(as.character(positions_with_read_ids$read_id_list[i]))
  #print(as.character(unique_resequenced_verified_reads$name_trimmed[tempindex]))
  #print('')
  if(unique_resequenced_verified_reads$paired[tempindex] == TRUE){ 
    positions_with_read_ids$paired[i]=TRUE
  }
  else { 
    positions_with_read_ids$paired[i]=FALSE
  }
}


all_positions_covered<-c()
all_positions_covered_readid<-c()

for(i in 1:length(positions_with_read_ids$positions_list[positions_with_read_ids$paired==TRUE])){ 
  print(i)
  if(positions_with_read_ids$paired[i]==TRUE){
  for(j in 1:36){
    temp<-(positions_with_read_ids$positions_list[i] + j)
    all_positions_covered<-append(all_positions_covered, temp)
    all_positions_covered_readid<-append(all_positions_covered_readid, as.character(positions_with_read_ids$read_id_list[i]))
    }
  }
}


all_positions_covered_df<-data.frame(all_positions_covered, all_positions_covered_readid)

ggplot(all_positions_covered_df) + 
  geom_freqpoly(aes(x = all_positions_covered_df$all_positions_covered), binwidth = 1) + 
  xlab('position') + 
  ylab('depth') + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,235000), breaks = c(0,47000,94000,141000,188000,235000)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1.5)) + 
  scale_y_continuous(expand = c(0, 0)) 
#xlim(0,250000) + 

write.csv(positions_with_read_ids,'positions_with_read_ids.csv')
