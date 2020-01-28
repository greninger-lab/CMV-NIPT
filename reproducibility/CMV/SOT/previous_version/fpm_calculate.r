#rebuild
#running after including full NT in blast db 

library(Rsamtools)
library(seqinr)
library(svMisc)
library(Biostrings)
library(svMisc)
library(ggplot2)
library(xlsx)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/')

blasthitsfile<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/deduplicated/blast_hits.csv')

#####edits
#blasthitsfile<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/deduplicated/blast_hits.csv')
####

#human_blasthitsfile<-read.csv('human_filtered_blast_hits.csv')
#blasthitsfile<-read.csv('cmv_blast_hits.csv')
#fastafile<-readDNAStringSet('/Users/gerbix/Documents/vikas/NIPT/21419_download/cmv_combined_masked.fasta')
#fastafile<-readDNAStringSet('/Users/gerbix/Documents/vikas/NIPT/new_samples/fastas/cmv_combined_masked.fasta')

####edits
fastafile<-readDNAStringSet('/Users/gerbix/Documents/vikas/NIPT/new_samples/deduplicated/cmv_combined_masked.fasta')
######

colnames(blasthitsfile)[1]<-'readname'
colnames(blasthitsfile)[2]<-'count'

read_counts_all<-read.csv2('/Users/gerbix/Documents/vikas/NIPT/new_samples/read_counts.csv', sep = ',',header = FALSE, col.names = c('sample', 'count'))
read_counts_all$sample<-as.character(read_counts_all$sample)
for(i in 1:length(read_counts_all$sample)){ 
  read_counts_all$sample[i]<-strsplit(read_counts_all$sample[i],'[.]')[[1]][1]
}
#####

blasthitsfile$readname<-as.character(blasthitsfile$readname)
blasthitsfile$unique_identifier<-NA
blasthitsfile$full_readname<-blasthitsfile$readname

blasthitsfile$sample_id<-NA
for(i in 1:nrow(blasthitsfile)){
  blasthitsfile$unique_identifier[i]<-strsplit(blasthitsfile$readname[i],'[-]')[[1]][2]
  blasthitsfile$readname[i]<-strsplit(blasthitsfile$readname[i],'[-]')[[1]][1]
  }
for(i in 1:nrow(blasthitsfile)){
  blasthitsfile$sample_id[i]<-strsplit(blasthitsfile$readname[i],'[.]')[[1]][1]
  }
blasthitsfile$fragments<-1
fragment_indexes<-which(duplicated(blasthitsfile$unique_identifier))
for(i in 1:length(fragment_indexes)){ 
  blasthitsfile$fragments[fragment_indexes[i]]<-0
  }



#####

#blank list for file in blast that have > 5 reads
positivefiles <- c() 
blasthitsfile$blast_pass<-FALSE
for( i in 1:nrow(blasthitsfile)){ 
  if(blasthitsfile$count[i] > 5) { 
    #print(blasthitsfile$readname[i])
    positivefiles<-append(positivefiles, as.character(blasthitsfile$full_readname[i]))
    blasthitsfile$blast_pass[i]<-TRUE
  }
}

#makes csv of all blast verified reads
nameslist<-c()
lengthslist<- c()
seqslist<-c()


for(i in 1:length(positivefiles)){ 
  progress(i, length(positivefiles))
  index=which(positivefiles[i]==names(fastafile))
  lengthslist = append(lengthslist, nchar(paste(fastafile[index])))
  nameslist<-append(nameslist,names(fastafile[index]))
  seqslist<-append(seqslist,paste(fastafile[index]))
  }
seqdata<-data.frame(nameslist,seqslist,lengthslist)
write.csv(seqdata,'blast_positive_sequence_info.csv')


positivestrimmed <- c()
for (word in positivefiles) { 
  word <- strsplit(word, '.sam')[[1]][1]
  positivestrimmed <- append(positivestrimmed, word )
}
unique_positivestrimmed<-unique(positivestrimmed)

######################
uniquefiles = c()
positvereads= c()
intermediatepositives = c()
intermediatepositives_reads = c()
positives<-c()
positivereads<-c()
weak_positives<-c()
weak_positives_reads<-c()


x<-blasthitsfile
unique(x$unique_identifier)
blasthitsfile<-blasthitsfile[which(!(duplicated(blasthitsfile$unique_identifier))),]

unique_file_list<-unique(blasthitsfile$sample_id)



for (i in 1:length(unique_file_list)) { 
  print(100*i/length(unique_file_list))
  count = 0
  identifier = unique_file_list[i]
  for ( j in 1:length(blasthitsfile$sample_id)) { 
    if(blasthitsfile$sample_id[j] == unique_file_list[i] && blasthitsfile$blast_pass[j]==TRUE) { 
      #print(i)
      #count = count + 1 
      count = count + blasthitsfile$fragments[j]
      #count = count + 1 + blasthitsfile$fragments[tempfragmentcount]
      #count = count + 1
      #print(blasthitsfile$fragments[tempfragmentcount])
      #print(tempfragmentcount)
      #print(identifier)
      tempidentifier<-blasthitsfile$sample_id[j]
    }
  }
  #print(identifier)
  # if( count > 15){ 
  #   print('whatthefuck')
  #   }
  if( count > 2) { 
    uniquefiles = append(uniquefiles,tempidentifier)
  }
  if(count < 15 && count >1) { 
    #print(identifier)
    intermediatepositives = append(intermediatepositives , tempidentifier)
    intermediatepositives_reads = append(intermediatepositives_reads , count)
  }
  if(count > 15) { 
    #print(identifier)
    positives = append(positives , tempidentifier)
    positivereads = append(positivereads , count)
  }
  if(count==1){ 
    weak_positives = append(weak_positives, tempidentifier)
    weak_positives_reads= append(weak_positives_reads, count)
  }
}





read_counts_all<-read_counts_all[complete.cases(read_counts_all$sample),]
#######################
#rpkm and rpm calculation for weak positives (1 read)
weakpositivesdf<-data.frame(weak_positives,weak_positives_reads)
weakpositivesdf <- weakpositivesdf[!duplicated(weakpositivesdf$weak_positives), ]
weakpositivesdf$rpkm<-NA
weakpositivesdf$read_counts<-NA
weakpositivesdf$rpm<-NA
for(i in 1:nrow(weakpositivesdf)){ 
  for( j in 1:nrow(read_counts_all)){ 
    if(read_counts_all$sample[j]==as.character(weakpositivesdf$weak_positives[i])){ 
      weakpositivesdf$rpkm[i]<-(weakpositivesdf$weak_positives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1000000))
      print((weakpositivesdf$weak_positives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1e9)))
      weakpositivesdf$read_counts[i]<-read_counts_all$count[j]
      weakpositivesdf$rpm[i]<-((weakpositivesdf$weak_positives_reads[i])*1e6)/(read_counts_all$count[j])

    }
  }
}

#rpkm and rpm calculation for intermediate positives (between 2-15 reads)
intermediatepositivesdf <- data.frame(intermediatepositives, intermediatepositives_reads)
intermediatepositivesdf <- intermediatepositivesdf[!duplicated(intermediatepositivesdf$intermediatepositives), ]
intermediatepositivesdf$rpkm<-NA
intermediatepositivesdf$read_counts<-NA
intermediatepositivesdf$rpm<-NA
read_counts_all$sample<-as.character(read_counts_all$sample)
for(i in 1:nrow(intermediatepositivesdf)){ 
  for( j in 1:nrow(read_counts_all)){ 
    if(read_counts_all$sample[j]==as.character(intermediatepositivesdf$intermediatepositives[i])){ 
      intermediatepositivesdf$rpkm[i]<-(intermediatepositivesdf$intermediatepositives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1000000))
      print((intermediatepositivesdf$intermediatepositives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1e9)))
      intermediatepositivesdf$read_counts[i]<-read_counts_all$count[j]
      intermediatepositivesdf$rpm[i]<-((intermediatepositivesdf$intermediatepositives_reads[i])*1e6)/(read_counts_all$count[j])
      
    }
  }
}




#rpkm and rpm calculation for strong positives (over 15 reads)
positivesdf<-data.frame(positives,positivereads)
positivesdf <- positivesdf[!duplicated(positivesdf$positives), ]
positivesdf$rpkm<-NA
positivesdf$read_counts<-NA
positivesdf$rpm<-NA
for(i in 1:nrow(positivesdf)){ 
  for( j in 1:nrow(read_counts_all)){ 
      if((strsplit(read_counts_all$sample[j],'[.]')[[1]][1])==as.character(positivesdf$positives[i])){ 
        positivesdf$rpkm[i]<-((positivesdf$positivereads[i])/((235646/1000)*(read_counts_all$count[j]/1000000)))
        print((positivesdf$positivereads[i])/((235646/1000)*(read_counts_all$count[j]/1e9)))
        positivesdf$read_counts[i]<-read_counts_all$count[j]
        positivesdf$rpm[i]<-((positivesdf$positivereads[i])*1e6)/(read_counts_all$count[j])
      }
    }
}

actuallyunique <- unique(uniquefiles)
write.csv(actuallyunique, 'other_cmv_positives.csv')
write.csv(weakpositivesdf,'cmv_weak_positives_table.csv')
write.csv(intermediatepositivesdf, 'cmv_intermediate_positives_table.csv')
write.csv(positivesdf, 'cmv_strong_positives_table.csv')
#cmv_all_positives_table<-rbind(positivesdf,intermediatepositivesdf,weakpositivesdf)





#graph to find a nice arbitrary rpkm cutoff
positivesdf$classification<-'strong positives'
intermediatepositivesdf$classification<-'intermediate positives'
weakpositivesdf$classification<-'weak positives'

graph_weakdf<-weakpositivesdf
colnames(graph_weakdf)[1]<-'sample'
colnames(graph_weakdf)[2]<-'count'

graph_positivesdf<-positivesdf
colnames(graph_positivesdf)[1]<-'sample'
colnames(graph_positivesdf)[2]<-'count'

graph_intermediatesdf<-intermediatepositivesdf
colnames(graph_intermediatesdf)[1]<-'sample'
colnames(graph_intermediatesdf)[2]<-'count'


df_for_graph<-rbind(graph_positivesdf,graph_intermediatesdf,graph_weakdf)
#df_for_graph<-mapply(c,positivesdf,intermediatepositivesdf)

for(i in 1:nrow(df_for_graph)){ 
    if(df_for_graph$rpm[i] > .3){ 
      df_for_graph$classification[i]<-'strong positive'
      }
    if(df_for_graph$rpm[i] < .3){ 
      df_for_graph$classification[i]<-'intermediate positive'
      }
    }

write.csv(df_for_graph,'all_sample_data.csv')



rpkm_plot<-ggplot(df_for_graph, aes(x=classification, y=rpkm)) + 
  geom_hline(yintercept = min(df_for_graph$rpkm[df_for_graph$classification=='strong positives'])) + 
  geom_hline(yintercept = max(df_for_graph$rpkm[df_for_graph$classification=='intermediate positives']), linetype = 'dashed') + 
  scale_y_continuous(trans = 'log10') +
  geom_point()
rpkm_plot
ggsave('CMV_logrpkm_by_classification.pdf', plot = last_plot())

reads_plot<-ggplot(df_for_graph, aes(x=classification, y=count)) + 
  geom_point() + 
  geom_hline(yintercept = min(df_for_graph$count[df_for_graph$classification=='strong positives'])) + 
  geom_hline(yintercept = max(df_for_graph$count[df_for_graph$classification=='intermediate positives']), linetype = 'dashed') + 
  scale_y_continuous(trans = 'log10') +
  ylab('# reads mapping to CMV')
reads_plot
ggsave('CMV_logcounts_by_classification.pdf', plot = last_plot())

write.xlsx(df_for_graph,'rpkm_values.xlsx')

rpkm_by_number_of_reads_plot<-ggplot(df_for_graph, aes(x=read_counts, y=rpkm, color=classification)) + 
  geom_hline(yintercept = min(df_for_graph$rpkm[df_for_graph$classification=='strong positives'])) + 
  geom_hline(yintercept = max(df_for_graph$rpkm[df_for_graph$classification=='intermediate positives']), linetype = 'dashed') + 
  scale_y_continuous(trans = 'log10') +
  xlab('number of reads in sample') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_point()
rpkm_by_number_of_reads_plot
ggsave('CMV_logrpkm_by_number_of_reads.pdf', plot = last_plot())




##same graphs but now using rpm


rpm_plot<-ggplot(df_for_graph, aes(x=classification, y=rpm)) + 
  geom_hline(yintercept = min(df_for_graph$rpm[df_for_graph$classification=='strong positives'])) + 
  geom_hline(yintercept = max(df_for_graph$rpm[df_for_graph$classification=='intermediate positives']), linetype = 'dashed') + 
  scale_y_continuous(trans = 'log10') +
  ylab('rpm of reads mapping to CMV') +
  geom_point()
rpm_plot
ggsave('CMV_logrpm_by_classification.pdf', plot = last_plot())

reads_plot<-ggplot(df_for_graph, aes(x=classification, y=count)) + 
  geom_point() + 
  geom_hline(yintercept = min(df_for_graph$count[df_for_graph$classification=='strong positives'])) + 
  geom_hline(yintercept = max(df_for_graph$count[df_for_graph$classification=='intermediate positives']), linetype = 'dashed') + 
  scale_y_continuous(trans = 'log10') +
  ylab('# reads mapping to CMV')
reads_plot
ggsave('CMV_logcounts_by_classification.pdf', plot = last_plot())

write.xlsx(df_for_graph,'rpkm_values.xlsx')

rpm_by_number_of_reads_plot<-ggplot(df_for_graph, aes(x=read_counts, y=rpm, color=classification)) + 
  geom_hline(yintercept = min(df_for_graph$rpm[df_for_graph$classification=='strong positives'])) + 
  geom_hline(yintercept = max(df_for_graph$rpm[df_for_graph$classification=='intermediate positives']), linetype = 'dashed') + 
  scale_y_continuous(trans = 'log10') +
  xlab('number of reads in sample') +
  ylab('rpm of reads mapping to CMV') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_point()
rpm_by_number_of_reads_plot
ggsave('CMV_logrpm_by_number_of_reads.pdf', plot = last_plot())











##filling in qpcr data for the all_data spreadsheet

edited_all_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/21419_download/all_sample_data.csv', header = TRUE)
edited_all_data$sample<-as.character(edited_all_data$sample)
qpcr_results<-read_excel("qpcr_results.xls", skip = 1)
qpcr_results$`Sample Name`<-gsub(pattern = '-',replacement = '',x=qpcr_results$`Sample Name`)

for(i in 1:nrow(qpcr_results)){ 
  for(j in 1:nrow(edited_all_data)){ 
    if(grepl(qpcr_results$`Sample Name`[i],edited_all_data$sample[j])){
      edited_all_data$cmv.quant[j]<-qpcr_results$`Quantity(10UL DNA)`[i]
      edited_all_data$bglobin.quant[j]<-qpcr_results$`Quantity(10UL DNA)__1`[i]
    }
    
    }
  }

write.csv(edited_all_data, 'all_data_combined.csv')





