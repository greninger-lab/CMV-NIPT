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

resequenced_verified_reads<-ReadFasta('/Users/gerbix/Documents/vikas/NIPT/new_samples/fastas/cmv_combined_masked.fasta')

#coverage graph: 

resequenced_bam<-scanBam('/Users/gerbix/Documents/vikas/NIPT/new_samples/merged_bam/new_cmv_merged.bam')

resequenced_verified_reads$name_trimmed<-NA
for(i in 1:nrow(resequenced_verified_reads)){ 
  progress(i,nrow(resequenced_verified_reads))
  resequenced_verified_reads$name_trimmed[i]<-strsplit(as.character(resequenced_verified_reads$name[i]),'-')[[1]][2]
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


resequenced_verified_reads$paired<-NA
for(i in 1:nrow(resequenced_verified_reads)){ 
  progress(i)
  tempindex<-which(grepl(as.character(resequenced_verified_reads$name_trimmed[i]),resequenced_bam[[1]]$qname))
  if(is.na(tempindex[2])){ 
    #print(index)
    resequenced_verified_reads$paired[i]<-FALSE
  }
  else{
    #print(tempindex)
    resequenced_verified_reads$paired[i]<-TRUE
  }
}





positions_list<-c()
read_id_list<-c()
for(i in 1:length(resequenced_verified_reads$name_trimmed)){
  index<- which(grepl(as.character(resequenced_verified_reads$name_trimmed[i]) , resequenced_bam[[1]]$qname))
  if(resequenced_verified_reads$paired[i]==FALSE){ 
     print((index))
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



# for(i in 1:nrow(resequenced_verified_reads)){ 
#   resequenced_verified_reads$pos[i]<-resequenced_bam[[1]]$pos[match_indices[i]]
#   }
# 
all_positions_covered<-c()
all_positions_covered_readid<-c()

for(i in 1:length(positions_with_read_ids$positions_list)){ 
  print(i)
  for(j in 1:36){
    temp<-(positions_with_read_ids$positions_list[i] + j)
    all_positions_covered<-append(all_positions_covered, temp)
    all_positions_covered_readid<-append(all_positions_covered_readid, as.character(positions_with_read_ids$read_id_list[i]))
  }
}




all_positions_covered_df<-data.frame(all_positions_covered, all_positions_covered_readid)

coverage_plot<-ggplot(all_positions_covered_df, aes(x = all_positions_covered_df$all_positions_covered)) + 
  geom_freqpoly( bins = 50) + 
  theme_classic() 
coverage_plot

ggsave('cmv_coverage_plot.pdf', coverage_plot,width = 3, height = 3)

coverage_freq_table<-data.frame(table(all_positions_covered_df$all_positions_covered))


positions_covered_by_freq<-ggplot(coverage_freq_table, aes(x = coverage_freq_table$Var1, y = coverage_freq_table$Freq, group = 1 )) + 
  geom_line() + 
  theme_classic() 

positions_covered_by_freq


coverage_plot<-ggplot(coverage_freq_table, aes(x = coverage_freq_table$Var1, y = coverage_freq_table$Freq, group = 1 )) + 
  geom_freqpoly(bins = 50, stat = 'identity') + 
  theme_classic() 
coverage_plot

#mapping = aes(x = x, y = y), stat = "identity"

x<-ggplot(all_positions_covered_df) + 
  geom_freqpoly(aes(x = all_positions_covered_df$all_positions_covered), binwidth = 1) + 
  xlab('position') + 
  ylab('depth') + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,235000), breaks = c(0,47000,94000,141000,188000,235000)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = .1)) + 
  scale_y_continuous(expand = c(0, 0)) 
  #xlim(0,250000) + 
x
ggsave(plot=x, "new_cmv_coverage.pdf", height = 3, width = 6)






