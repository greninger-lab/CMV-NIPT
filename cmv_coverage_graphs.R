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
#to_remove<- which(duplicated(as.character(resequenced_verified_reads$name_trimmed)))
#unique_resequenced_verified_reads<-resequenced_verified_reads[-to_remove,]

resequenced_bam_names<-resequenced_bam[[1]]$qname
# resequenced_bam_names_trimmed<-c()
# for(i in 1:length(resequenced_bam_names)){ 
#   progress(i,nrow(verified_reads))
#   resequenced_bam_names_trimmed <-append(resequenced_bam_names_trimmed,strsplit(as.character(resequenced_bam_names[i]),'/')[[1]][1])
# }
# 

match_indices<-c()
for(i in 1:nrow(resequenced_verified_reads)){ 
  match_indices<-append(match_indices, which(as.character(resequenced_verified_reads$name_trimmed[i]) == resequenced_bam_names))
  }


unique_resequenced_verified_reads$pos<-NA

for(i in 1:nrow(unique_resequenced_verified_reads)){ 
  unique_resequenced_verified_reads$pos[i]<-resequenced_bam[[1]]$pos[match_indices[i]]
  }

all_positions_covered<-c()
all_positions_covered_readid<-c()

for(i in 1:length(unique_resequenced_verified_reads$pos)){ 
  print(i)
  for(j in 1:36){
    temp<-(unique_resequenced_verified_reads$pos[i] + j)
  all_positions_covered<-append(all_positions_covered, temp)
  all_positions_covered_readid<-append(all_positions_covered_readid, as.character(unique_resequenced_verified_reads$name_trimmed[i]))
  
    }
  }

all_positions_covered_df<-data.frame(all_positions_covered, all_positions_covered_readid)

coverage_plot<-ggplot(all_positions_covered_df, aes(x = all_positions_covered_df$all_positions_covered)) + 
  geom_freqpoly(bins = 70) + 
  theme_classic() 
coverage_plot

ggsave('cmv_coverage_plot.pdf', coverage_plot,width = 3, height = 3)







