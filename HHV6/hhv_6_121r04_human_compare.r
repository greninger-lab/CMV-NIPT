library(ggplot2)

###121r04 human 
setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/figure_2_final/71719')
#read_counts <- read.csv("~/Documents/vikas/NIPT/clip_removed/cmv_full/read_counts.csv")
filenames = '/Volumes/Seagate8Tb1/nipt/bams/bams_for_vikas/3P13_E02_CFFv1_NB0015.final.bam.read_length_counts.txt'
plotslist<-c()


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


human_frequencies$isize<-as.numeric(as.character(human_frequencies$isize))
x<-human_frequencies[order(human_frequencies$isize),]

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


human_plot<-ggplot(human_frequencies, aes( x = human_frequencies$isize , y = human_frequencies$percent)) + 
  geom_vline(xintercept = human_median) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
human_plot
ggsave(plot = human_plot, 'Resequenced_human_fragment_length.pdf')


###hhv6 isizes
hhv_6_combined<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/figure_2_final/71719/hhv6_isizes_all.csv')



plot<-ggplot(hhv_6_combined, aes( x=hhv_6_combined$length, y= hhv_6_combined$percent)) + 
  geom_histogram(binwidth = 5, aes(fill = combineddf$sample, y=(100 * ..count../sum(..count..)))) + 
  geom_line(data = , aes(x = A, y = B), col = "red") +
  xlim(0,500) + 
  ylab('percent')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")
  
plot



hhv6_ks_test<-ks.test(human_isize_expanded,hhv_6_combined$lengthlist)
hhv6_ks_test


hhv_6_isize<-combined$length
hhv_6_isize$type<-'hhv6'

human_isize<-data.frame(human_isize_expanded)
colnames(human_isize)[1]<-'length'
human_isize$type<-'Human'

human_hhv6_combined<-rbind(human_isize,hhv_6_isize)



cum_frequency<-ggplot(subsampled_df, aes(x = subsampled_df$isizes, color = subsampled_df$type)) + 
  theme_classic() +  
  theme(legend.position='none') + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') + 
  stat_ecdf(geom = 'step', size  =1 ) 
cum_frequency


plot(ecdf(combined$length))
plot(ecdf(human_isize_expanded), add = TRUE, col = 2 )
plot(ecdf(cmv))





uniques<-unique(combineddf$sample)

for(i in 1:length(uniques)){ 
  print(uniques[i])
  print( nrow(combineddf[combineddf$sample==uniques[i],]))
  }

#99P03_C01
#120R22_F03
#104P04_D01
#98p11
#97P10_B02
#93R20_D03

to_remove<-c('99P03_C01','120R22_F03','104P04_D01','98P11','97P10_B02','93R20_D03')


lowclusterremoved<-combineddf[-which(grepl(paste(to_remove,collapse="|"),combineddf$sample)),]

plot(ecdf(lowclusterremoved$length))
plot(ecdf(human_isize_expanded), add = TRUE, col = 2 )


lowclusteronly<-combineddf[which(grepl(paste(to_remove,collapse="|"),combineddf$sample)),]

plot(ecdf(lowclusterremoved$length))
plot(ecdf(human_isize_expanded), add = TRUE, col = 2 )
plot(ecdf(lowclusteronly$length), add = TRUE, col = 3 )


ks.test(human_isize_expanded, combineddf$length)


