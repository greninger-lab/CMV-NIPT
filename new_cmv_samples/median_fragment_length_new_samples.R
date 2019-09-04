library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)

setwd('/Users/gerbix/Documents/vikas/NIPT/new_samples')
#read_counts <- read.csv("~/Documents/vikas/NIPT/clip_removed/cmv_full/read_counts.csv")
blast_hits_file<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/blast_hits.csv', header = FALSE, col.names = c('full_ID','count'))
filenames = list.files('/Users/gerbix/Documents/vikas/NIPT/new_samples',pattern = '.bam$')
plotslist<-c()


blast_hits_file$full_ID<-as.character(blast_hits_file$full_ID)
for( i in 1:nrow(blast_hits_file)){
  blast_hits_file$read_ID[i]<-strsplit(blast_hits_file$full_ID, '-')[[i]][2]
} 

blast_hits_file<-blast_hits_file[blast_hits_file$count>75,]
blast_hits_file<-blast_hits_file[-(which(duplicated(blast_hits_file$read_ID))),]

isize<-c()
read_id<-c()
sample<-c()
for(i in 1:length(filenames)){ 
  #print(i)
  temp_bam<-scanBam(filenames[i])
  if((identical(temp_bam[[1]]$qname, character(0)))){ 
    print(i)
    next
  }
  all_bam_read_IDs<-temp_bam[[1]]$qname  
  matches<-c()
  for(j in 1:nrow(blast_hits_file)){ 
    #progress(j)
    for(k in 1:length(all_bam_read_IDs)){
      print(blast_hits_file$read_ID[j])
      if( !(is.na(blast_hits_file$read_ID[j])) & grepl( blast_hits_file$read_ID[j], all_bam_read_IDs[k])) { 
        isize<-append(isize, temp_bam[[1]]$isize[k])
        read_id<-append(read_id, paste(temp_bam[[1]]$qname[k],'.1'))
        sample<-append(sample, filenames[i])
      }
    }
  }
}
combined<-data.frame(sample, read_id, isize)
combined<-combined[combined$isize>0,]
combined<-combined[complete.cases(combined$isize),]
cmv_plot<-ggplot(combined, aes( x=combined$isize)) + 
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(combined$isize)) +
  xlim(0,500) + 
  #ylim(0,200) + 
  ylab('occurences')+
  xlab('fragment length')+
  #annotate("text", x = 400, y = 18 , label =  paste0('mean=', mean(combined$isize))) + 
  annotate("text", x = 400, y = 15, label =  paste0('median=', median(combined$isize)))  +
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")
cmv_plot


for(i in (unique(combined$sample))){ 
  print(i)
  print(median(combined$isize[combined$sample==i & combined$isize < 500]))
  print('')
    }


P17<-combined[combined$sample=='244P17_A03_CFFv2_NB0289.sam.bam' & combined$isize < 500,]
median(P17$isize[P17$isize<500])






















setwd('/Users/gerbix/Documents/vikas/NIPT/new_samples')
read_counts <- read.csv("/Users/gerbix/Documents/vikas/NIPT/new_samples/read_counts.csv")
filenames = list.files(pattern = '*.bam.results.txt$')
plotslist<-c()

savePlot <- function(myPlot, plotname) {
  pdf(plotname, width = 8.5, height = 11)
  print(myPlot)
  dev.off()
}

dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(filenames)){ 
  print(i)
  #print(filenames[i])
  if (file.size(filenames[i]) != 0) {
    histo<- read.table(filenames[i], quote="\"", comment.char="")
    colnames(histo)[2]<-'length'
    colnames(histo)[1]<-'occurences'
    histo$type<-strsplit(filenames[i],".sam")[[1]][1]
    histo<-histo[histo$length>0,]
    if(nrow(histo)>0){
      toremove<-which(histo[,2]> 900)
      if(length(toremove)>0){ 
        histo<-histo[-toremove,]
      }
      totallengths= append(totallengths,histo$length)
      totaloccurences =append(totaloccurences,histo$occurences)
      dflist[[i]]<-histo
    }
  }
}

combineddf <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()
read_counts$sample<-as.character(read_counts$sample)

count = 0 
for(i in 1:nrow(combineddf)){ 
  temp = combineddf$occurences[i] * combineddf$length[i]
  count = count + temp
}
human_mean = count/ sum(combineddf$occurences)

median_list<-c()
for(i in 1:nrow(combineddf)){ 
  progress(i, nrow(combineddf))
  median_list<-append(median_list,rep(combineddf$length[i],combineddf$occurences[i]))
}
human_median<-median(median_list)

human_plot<-ggplot(combineddf, aes( x=combineddf$length, y = combineddf$occurences, color = combineddf$type)) + 
  geom_point() + 
  #geom_histogram(binwidth = 5) + 
  #  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = human_median) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="none")  +
  annotate("text", x = 400, y = 1500000 , label =  paste0('mean=', mean(human_mean))) +
  annotate("text", x = 400, y = 1300000, label =  paste0('median=', median(human_median)))  + 
  ggtitle('Human median fragment length')

human_plot


