library(ggplot2)

#pulled human insert sizes from NIPT CMV sample 3P13
setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6')
filenames = '/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/3P13_read_lengths.txt'

dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(filenames)){ 
  print(i)
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

human_frequencies$percent<-(100 * human_frequencies$freq) / (sum(human_frequencies$freq))

#human median calculation
human_isize_expanded<-rep(human_frequencies$isize, human_frequencies$freq) 
human_isize_expanded_below250<-human_isize_expanded[human_isize_expanded < 500]
human_median<-median(as.numeric(as.character(human_isize_expanded_below250)))


human_frequencies$isize<-as.numeric(as.character(human_frequencies$isize))


###hhv6 isizes
hhv_6_combined<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/hhv6_isizes_all.csv')
 
hhv6_ks_test<-ks.test(human_isize_expanded,hhv_6_combined$lengthlist)
hhv6_ks_test


hhv_6_isize<-data.frame(hhv_6_combined$lengthlist)
colnames(hhv_6_isize)[1]<-'length'
hhv_6_isize$type<-'hhv6'

human_isize<-data.frame(human_isize_expanded)
colnames(human_isize)[1]<-'length'
human_isize$type<-'Human'

human_hhv6_combined<-rbind(human_isize,hhv_6_isize)


cum_frequency<-ggplot(human_hhv6_combined, aes(x = human_hhv6_combined$length, color = human_hhv6_combined$type)) + 
  theme_classic() +  
  theme(legend.position='none') + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') + 
  stat_ecdf(geom = 'step', size  =1 ) 
cum_frequency


plot(ecdf(combined$length))
plot(ecdf(human_isize_expanded), add = TRUE, col = 2 )
plot(ecdf(cmv))

for(i in 1:nrow(human_hhv6_combined)){ 
  if(human_hhv6_combined$type=='H')
}

uniques<-unique(combineddf$sample)

for(i in 1:length(uniques)){ 
  print(uniques[i])
  print( nrow(combineddf[combineddf$sample==uniques[i],]))
}

to_remove<-c('99P03_C01','120R22_F03','104P04_D01','98P11','97P10_B02','93R20_D03')

lowclusterremoved<-combineddf[-which(grepl(paste(to_remove,collapse="|"),combineddf$sample)),]

plot(ecdf(lowclusterremoved$length))
plot(ecdf(human_isize_expanded), add = TRUE, col = 2 )


lowclusteronly<-combineddf[which(grepl(paste(to_remove,collapse="|"),combineddf$sample)),]

plot(ecdf(lowclusterremoved$length))
plot(ecdf(human_isize_expanded), add = TRUE, col = 2 )
plot(ecdf(lowclusteronly$length), add = TRUE, col = 3 )


ks.test(human_isize_expanded, combineddf$length)

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

# saving plot crashes R 
# ggsave(cum_frequency, 'figure_5C.pdf', height = 3 , width = 3)


