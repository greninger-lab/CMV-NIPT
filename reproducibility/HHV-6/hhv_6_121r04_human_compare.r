library(ggplot2)

###121r04 human 
setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6')
dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(filenames)){ 
  print(i)
  #print(filenames[i])
  if (file.size(filenames[i]) != 0) {
    histo<- read.table(filenames[i], quote="\"", comment.char="")
    colnames(histo)[2]<-'isize'
    colnames(histo)[1]<-'freq'
    histo$type<-'Human'
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

human_frequencies <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()
read_counts$sample<-as.character(read_counts$sample)

human_frequencies$percent<-(100 * human_frequencies$freq) / (sum(human_frequencies$freq))

#human median calculation
human_isize_expanded<-rep(human_frequencies$isize, human_frequencies$freq) 
human_isize_expanded_below250<-human_isize_expanded[human_isize_expanded < 500]
human_median<-median(as.numeric(as.character(human_isize_expanded_below250)))


hhv6_isizes_all<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/hhv6_isizes_all.csv')

human_isizes_all<-data.frame(human_isize$length)
colnames(human_isizes_all)[1]<-'lengthlist'
human_isizes_all$X<-NA
human_isizes_all$samplelist<-'Human'
human_isizes_all$sample_trimmed<-'Human'
human_isizes_all$cluster<-'Human'

all_combined<-rbind(human_isizes_all, hhv6_isizes_all)

cum_frequency<-ggplot(all_combined, aes(x = all_combined$lengthlist, color = all_combined$cluster)) + 
  stat_ecdf(geom = 'step', size  =1 , show.legend = TRUE) + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlim(0,500) + 
  theme(legend.position='bottom') + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') 
cum_frequency
ggsave(cum_frequency, 'figure_5C.pdf', height = 3 , width = 3)


