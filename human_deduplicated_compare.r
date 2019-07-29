library(ggplot2)
#comparing human with and without duplicates 


setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating')
deduplicated_human = '/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/121r04_original_human_duplicates_removed_read_lengths_count.txt'
plotslist<-c()


dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(deduplicated_human)){ 
  print(i)
  #print(filenames[i])
  if (file.size(deduplicated_human[i]) != 0) {
    histo<- read.table(deduplicated_human[i], quote="\"", comment.char="")
    colnames(histo)[2]<-'isize'
    colnames(histo)[1]<-'freq'
    histo$type<-'deduplicated'
    histo<-histo[histo$isize>0,]
    if(nrow(histo)>0){
      toremove<-which(histo[,2]> 500)
      if(length(toremove)>0){ 
        histo<-histo[-toremove,]
      }
      totallengths= append(totallengths,histo$isize)
      totaloccurences =append(totaloccurences,histo$freq)
      dflist[[i]]<-histo
    }
  }
}

deduplicated_human_frequencies <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()

deduplicated_human_frequencies$percent<-(100 * deduplicated_human_frequencies$freq) / (sum(deduplicated_human_frequencies$freq))

#human median calculation
human_isize_expanded_deduplicated<-rep(deduplicated_human_frequencies$isize, deduplicated_human_frequencies$freq) 
human_isize_expanded_below500<-human_isize_expanded_deduplicated[human_isize_expanded_deduplicated < 500]
human_median_deduplicated<-median(as.numeric(as.character(human_isize_expanded_below500)))


deduplicated_human_frequencies$isize<-as.numeric(as.character(deduplicated_human_frequencies$isize))
x<-deduplicated_human_frequencies[order(deduplicated_human_frequencies$isize),]

x<-x[x$isize < 250 ,]

count = 1 
sum = 0 
x$percent<- (100 * x$freq)/ sum(x$freq)
while( sum < 50) { 
  print(sum)
  print(x$isize[count])
  sum = sum +  x$percent[count]
  count = count + 1 
}


human_plot<-ggplot(deduplicated_human_frequencies, aes( x = deduplicated_human_frequencies$isize , y = deduplicated_human_frequencies$percent)) + 
  geom_vline(xintercept = human_median) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
human_plot
ggsave(plot = human_plot, 'Resequenced_human_fragment_length.pdf')









#Human with duplicates

setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating')
human_with_duplicates = '/Users/gerbix/Documents/vikas/NIPT/31119_download/figure_2_final/121R04_D01_CFFv1_NA0144.final.bam.results.txt'
plotslist<-c()


dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(human_with_duplicates)){ 
  print(i)
  #print(filenames[i])
  if (file.size(human_with_duplicates[i]) != 0) {
    histo<- read.table(human_with_duplicates[i], quote="\"", comment.char="")
    colnames(histo)[2]<-'isize'
    colnames(histo)[1]<-'freq'
    histo$type<-'with_duplicates'
    histo<-histo[histo$isize>0,]
    if(nrow(histo)>0){
      toremove<-which(histo[,2]> 500)
      if(length(toremove)>0){ 
        histo<-histo[-toremove,]
      }
      totallengths= append(totallengths,histo$isize)
      totaloccurences =append(totaloccurences,histo$freq)
      dflist[[i]]<-histo
    }
  }
}

human_with_duplicates <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()

human_with_duplicates$percent<-(100 * human_with_duplicates$freq) / (sum(human_with_duplicates$freq))

#human median calculation
human_isize_expanded_duplicates<-rep(human_with_duplicates$isize, human_with_duplicates$freq) 
human_isize_expanded_below500<-human_isize_expanded_duplicates[human_isize_expanded_duplicates < 500]
human_median_duplicates<-median(as.numeric(as.character(human_isize_expanded_below500)))


human_with_duplicates$isize<-as.numeric(as.character(human_with_duplicates$isize))
x<-human_with_duplicates[order(human_with_duplicates$isize),]

x<-x[x$isize < 250 ,]

count = 1 
sum = 0 
x$percent<- (100 * x$freq)/ sum(x$freq)
while( sum < 50) { 
  print(sum)
  print(x$isize[count])
  sum = sum +  x$percent[count]
  count = count + 1 
}


human_plot<-ggplot(human_with_duplicates, aes( x = human_with_duplicates$isize , y = human_with_duplicates$percent)) + 
  geom_vline(xintercept = human_median) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
human_plot
ggsave(plot = human_plot, 'Resequenced_human_fragment_length.pdf')



####only the removed duplicates
extracted_duplicates<-'/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/duplicates_read_lengths.txt'
plotslist<-c()
dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(extracted_duplicates)){ 
  print(i)
  #print(filenames[i])
  if (file.size(extracted_duplicates[i]) != 0) {
    histo<- read.table(extracted_duplicates[i], quote="\"", comment.char="")
    colnames(histo)[2]<-'isize'
    colnames(histo)[1]<-'freq'
    histo$type<-'extracted_duplicates'
    histo<-histo[histo$isize>0,]
    if(nrow(histo)>0){
      toremove<-which(histo[,2]> 500)
      if(length(toremove)>0){ 
        histo<-histo[-toremove,]
      }
      totallengths= append(totallengths,histo$isize)
      totaloccurences =append(totaloccurences,histo$freq)
      dflist[[i]]<-histo
    }
  }
}

extracted_duplicates <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()

extracted_duplicates$percent<-(100 * extracted_duplicates$freq) / (sum(extracted_duplicates$freq))

#human median calculation
extractedduplicates_expanded<-rep(extracted_duplicates$isize, extracted_duplicates$freq) 
extractedduplicates_expanded_below500<-extractedduplicates_expanded[extractedduplicates_expanded < 500]
extractedduplicates_expanded_median<-median(as.numeric(as.character(extractedduplicates_expanded_below500)))


extracted_duplicates$isize<-as.numeric(as.character(extracted_duplicates$isize))
x<-extracted_duplicates[order(extracted_duplicates$isize),]

x<-x[x$isize < 250 ,]

count = 1 
sum = 0 
x$percent<- (100 * x$freq)/ sum(x$freq)
while( sum < 50) { 
  print(sum)
  print(x$isize[count])
  sum = sum +  x$percent[count]
  count = count + 1 
}


human_plot<-ggplot(extracted_duplicates, aes( x = extracted_duplicates$isize , y = extracted_duplicates$percent)) + 
  geom_vline(xintercept = human_median) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
human_plot











#combining duplicated and deduplicated
#human_with_duplicates$type<-'with duplicates'
#deduplicated_human_frequencie$type<-'without duplicates'

duplicated_deduplicated_combined<-rbind(deduplicated_human_frequencies, human_with_duplicates, extracted_duplicates)

duplication_distribution_plot<-ggplot(duplicated_deduplicated_combined, aes( x = duplicated_deduplicated_combined$isize , y = duplicated_deduplicated_combined$percent, color = duplicated_deduplicated_combined$type)) + 
  geom_vline(xintercept = 168) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
duplication_distribution_plot
ggsave(plot = human_plot, 'Resequenced_human_fragment_length.pdf')



duplicated_expaned<-data.frame(human_isize_expanded_duplicates)
duplicated_expaned$type<-'duplicated'
colnames(duplicated_expaned)<-c('isizes', 'type')



deduplicated_expanded<-data.frame(human_isize_expanded_deduplicated)
deduplicated_expanded$type<-'deduplicated'
colnames(deduplicated_expanded)<-c('isizes', 'type')


extracted_duplicates_expanded_df<-data.frame(extractedduplicates_expanded)
extracted_duplicates_expanded_df$type<-'extracted_duplicates'
colnames(extracted_duplicates_expanded_df)<-c('isizes', 'type')



expanded_combined<-rbind(duplicated_expaned, deduplicated_expanded,(extracted_duplicates_expanded_df))

duplication_distribution_histogram<-ggplot(expanded_combined, aes(x=expanded_combined$isizes, color =expanded_combined$type)) +
                                             geom_freqpoly(binwidth = 5 ) + 
  theme_classic() + 
  xlab('Insert sizes')

duplication_distribution_histogram


