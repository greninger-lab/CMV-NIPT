library(Rsamtools)
library(seqinr)
library(svMisc)
library(Biostrings)
library(svMisc)
library(ggplot2)
library(xlsx)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility')

#Provide the blast hit count table here 
blasthitsfile<-read.csv('/Users/gerbix/Documents/vikas/NIPT/all_deduplicated/blast_hits.csv')

#provide the repeatmaskedmasked fasta here
fastafile<-readDNAStringSet('/Users/gerbix/Documents/vikas/NIPT/all_deduplicated/original_bams_deduplicated/cmv_combined_masked.fasta')
colnames(blasthitsfile)[1]<-'readname'
colnames(blasthitsfile)[2]<-'count'

#provide a csv of all read counts in the original FASTQ files here (column1=name, column2=read_count)
read_counts_all<-read.csv2('/Users/gerbix/Documents/vikas/NIPT/31119_download/read_counts_all_deduplicated.csv', sep = ',')

#Creating a lookup table for later on 
read_counts_all$sample<-as.character(read_counts_all$sample)
for(i in 1:length(read_counts_all$sample)){ 
  read_counts_all$sample[i]<-strsplit(read_counts_all$sample[i],'[.]')[[1]][1]
}

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

uniquefiles = c()
positvereads= c()
intermediatepositives = c()
intermediatepositives_reads = c()
positives<-c()
positivereads<-c()
weak_positives<-c()
weak_positives_reads<-c()

#convert reads to fragments based on unique identifiers in the blast file
blasthitsfile<-blasthitsfile[which(!(duplicated(blasthitsfile$unique_identifier))),]
blasthitsfile<-blasthitsfile[-which(blasthitsfile$full_readname=='121R27_C04_CFFv1_NA0147.37.fasta-:23112:19889:5382'),]





#matches read IDs back to their original SAM files 
unique_file_list<-unique(blasthitsfile$sample_id)
for (i in 1:length(unique_file_list)) { 
  print(100*i/length(unique_file_list))
  count = 0
  identifier = unique_file_list[i]
  for ( j in 1:length(blasthitsfile$sample_id)) { 
    if(blasthitsfile$sample_id[j] == unique_file_list[i] && blasthitsfile$blast_pass[j]==TRUE) { 
      count = count + blasthitsfile$fragments[j]
      tempidentifier<-blasthitsfile$sample_id[j]
    }
  }
  if( count > 2) { 
    uniquefiles = append(uniquefiles,tempidentifier)
  }
  if(count < 15 && count >1) { 
    intermediatepositives = append(intermediatepositives , tempidentifier)
    intermediatepositives_reads = append(intermediatepositives_reads , count)
  }
  if(count > 15) { 
    positives = append(positives , tempidentifier)
    positivereads = append(positivereads , count)
  }
  if(count==1){ 
    weak_positives = append(weak_positives, tempidentifier)
    weak_positives_reads= append(weak_positives_reads, count)
  }
}






#fpkm and fpm calculation for weak positives (1 read)
weakpositivesdf<-data.frame(weak_positives,weak_positives_reads)
weakpositivesdf <- weakpositivesdf[!duplicated(weakpositivesdf$weak_positives), ]
weakpositivesdf$fpkm<-NA
weakpositivesdf$read_counts<-NA
weakpositivesdf$fpm<-NA
for(i in 1:nrow(weakpositivesdf)){ 
  for( j in 1:nrow(read_counts_all)){ 
    if(read_counts_all$sample[j]==as.character(weakpositivesdf$weak_positives[i])){ 
      weakpositivesdf$fpkm[i]<-(weakpositivesdf$weak_positives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1000000))
      print((weakpositivesdf$weak_positives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1e9)))
      weakpositivesdf$read_counts[i]<-read_counts_all$count[j]
      weakpositivesdf$fpm[i]<-((weakpositivesdf$weak_positives_reads[i])*1e6)/(read_counts_all$count[j])

    }
  }
}

#fpkm and fpm calculation for intermediate positives (between 2-15 reads)
intermediatepositivesdf <- data.frame(intermediatepositives, intermediatepositives_reads)
intermediatepositivesdf <- intermediatepositivesdf[!duplicated(intermediatepositivesdf$intermediatepositives), ]
intermediatepositivesdf$fpkm<-NA
intermediatepositivesdf$read_counts<-NA
intermediatepositivesdf$fpm<-NA
read_counts_all$sample<-as.character(read_counts_all$sample)
for(i in 1:nrow(intermediatepositivesdf)){ 
  for( j in 1:nrow(read_counts_all)){ 
    if(read_counts_all$sample[j]==as.character(intermediatepositivesdf$intermediatepositives[i])){ 
      intermediatepositivesdf$fpkm[i]<-(intermediatepositivesdf$intermediatepositives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1000000))
      print((intermediatepositivesdf$intermediatepositives_reads[i])/((235646/1000)*(read_counts_all$count[j]/1e9)))
      intermediatepositivesdf$read_counts[i]<-read_counts_all$count[j]
      intermediatepositivesdf$fpm[i]<-((intermediatepositivesdf$intermediatepositives_reads[i])*1e6)/(read_counts_all$count[j])
      
    }
  }
}




#fpkm and fpm calculation for strong positives (over 15 reads)
positivesdf<-data.frame(positives,positivereads)
positivesdf <- positivesdf[!duplicated(positivesdf$positives), ]
positivesdf$fpkm<-NA
positivesdf$read_counts<-NA
positivesdf$fpm<-NA
for(i in 1:nrow(positivesdf)){ 
  for( j in 1:nrow(read_counts_all)){ 
      if((strsplit(read_counts_all$sample[j],'[.]')[[1]][1])==as.character(positivesdf$positives[i])){ 
        positivesdf$fpkm[i]<-((positivesdf$positivereads[i])/((235646/1000)*(read_counts_all$count[j]/1000000)))
        print((positivesdf$positivereads[i])/((235646/1000)*(read_counts_all$count[j]/1e9)))
        positivesdf$read_counts[i]<-read_counts_all$count[j]
        positivesdf$fpm[i]<-((positivesdf$positivereads[i])*1e6)/(read_counts_all$count[j])
      }
    }
}

actuallyunique <- unique(uniquefiles)
cmv_all_positives_table<-rbind(positivesdf,intermediatepositivesdf,weakpositivesdf)





#graph to find a nice arbitrary fpkm cutoff
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


for(i in 1:nrow(df_for_graph)){ 
    if(df_for_graph$fpm[i] > .3){ 
      df_for_graph$classification[i]<-'strong positive'
      }
    if(df_for_graph$fpm[i] < .3){ 
      df_for_graph$classification[i]<-'intermediate positive'
      }
    }

write.csv(df_for_graph,'all_sample_data.csv')






