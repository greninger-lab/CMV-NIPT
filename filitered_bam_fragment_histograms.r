library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)

setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/filtered_bams')
#read_counts <- read.csv("~/Documents/vikas/NIPT/clip_removed/cmv_full/read_counts.csv")
blast_hits_file<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/blast_hits.csv', header = FALSE, col.names = c('full_ID','count'))
filenames = list.files('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/original_bams',pattern = '.bam$')
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
    if(grepl( blast_hits_file$read_ID[j], all_bam_read_IDs[k])) { 
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
  annotate("text", x = 400, y = 18 , label =  paste0('mean=', mean(combined$isize))) + 
  annotate("text", x = 400, y = 15, label =  paste0('median=', median(combined$isize)))  +
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")
cmv_plot


ggsave('combined_cmv_plot.pdf',cmv_plot)





#####human stuff######


setwd('/Volumes/Seagate8Tb1/nipt/bams/bams_realigned_to_hg38')
read_counts <- read.csv("~/Documents/vikas/NIPT/clip_removed/cmv_full/read_counts.csv")
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
setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches')

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

ggsave('human_plot.pdf', human_plot)




#overlaying plots

human_plot<-ggplot(combineddf, aes( x=combineddf$length, y = combineddf$occurences, color = combineddf$type)) +   
  geom_point(aes = ) + 
  #geom_histogram(binwidth = 5) + 
  #  geom_histogram(binwidth = 5) +
  #geom_vline(xintercept = median(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="none") 
#annotate("text", x = 400, y = max(histo$occurences)-5, label =  paste0('mean=', mean(combineddf$length))) + 
#annotate("text", x = 400, y = max(histo$occurences)-7, label =  paste0('median=', median(combineddf$length)))
human_plot







###121R04 original plot
R04_original<-combined[combined$sample=='121R04_D01_CFFv1_NA0144.sam.bam',]
r04_original_plot<-ggplot(R04_original, aes( x=R04_original$isize)) + 
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(R04_original$isize)) +
  xlim(0,500) + 
  #ylim(0,200) + 
  ylab('occurences')+
  xlab('fragment length')+
  annotate("text", x = 400, y = 18 , label =  paste0('mean=', mean(R04_original$isize))) + 
  annotate("text", x = 400, y = 15, label =  paste0('median=', median(R04_original$isize)))  +
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom") + 
  ggtitle('121R04 original')
r04_original_plot
ggsave('121R04_original.pdf', r04_original_plot)






###12104 resequenced 
setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/resequenced_sample')
blast_hits_file<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/resequenced_sample/resequenced_cmv_vs_full_nt.txt_blast_hits.csv', header = FALSE, col.names = c('full_ID','count','x'))
r04_resequenced_BamFile<-scanBam('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/resequenced_sample/resequenced_local_cmv_aligned.bam')
duplicated_file_data<-read.csv('positions_with_read_ids.csv')
plotslist<-c()


blast_hits_file$full_ID<-as.character(blast_hits_file$full_ID)
for( i in 1:nrow(blast_hits_file)){
  progress(i,nrow(blast_hits_file))
  blast_hits_file$read_ID[i]<-strsplit(blast_hits_file$full_ID, '-')[[i]][2]
} 

#blast_hits_file<-blast_hits_file[blast_hits_file$count>75,]
blast_hits_file<-blast_hits_file[-(which(duplicated(blast_hits_file$read_ID))),]

isize<-c()
read_id<-c()
sample<-c()
  all_bam_read_IDs<-r04_resequenced_BamFile[[1]]$qname  
  matches<-c()
  for(j in 1:nrow(blast_hits_file)){ 
    progress(j,nrow(blast_hits_file))
    for(k in 1:length(all_bam_read_IDs)){
      #print(grepl( blast_hits_file$read_ID[j], all_bam_read_IDs[k])) 
        
      if(grepl( blast_hits_file$read_ID[j], all_bam_read_IDs[k])) { 
       # print ('ok')
        isize<-append(isize, r04_resequenced_BamFile[[1]]$isize[k])
        #print(isize)
        read_id<-append(read_id, paste(r04_resequenced_BamFile[[1]]$qname[k],'.1'))
        sample<-append(sample, '121r04_resequenced')
      }
    }
  }
r04_resequenced_combined<-data.frame(sample, read_id, isize)
r04_resequenced_combined<-r04_resequenced_combined[r04_resequenced_combined$isize>0,]
r04_resequenced_combined<-r04_resequenced_combined[complete.cases(r04_resequenced_combined$isize),]
r04_resequenced_plot<-ggplot(r04_resequenced_combined, aes( x=combined$isize)) + 
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(r04_resequenced_combined$isize)) +
  xlim(0,500) + 
  #ylim(0,200) + 
  ylab('occurences')+
  xlab('fragment length')+
  annotate("text", x = 400, y = 100 , label =  paste0('mean=', mean(r04_resequenced_combined$isize))) + 
  annotate("text", x = 400, y = 90, label =  paste0('median=', median(r04_resequenced_combined$isize)))  +
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom") + 
  ggtitle('121R04 resequenced')
r04_resequenced_plot


ggsave('121r04_resequenced_plot.pdf',r04_resequenced_plot)


write.csv(combined,'121r04_resequenced.csv')


not_unique<-which(duplicated(duplicated_file_data$positions_list))
duplicated_file_data<-duplicated_file_data[-not_unique,]
#duplicated_file_data<-duplicated_file_data[complete.cases(duplicated_file_data$positions_list),]

duplicates_to_keep<-c()
#non_paired<-duplicated_file_data$read_id_list[duplicated_file_data$paired==TRUE]
for(i in 1:nrow(duplicated_file_data)){ 
  if(duplicated_file_data$paired[i]==TRUE){ 
  for(j in 1:nrow(r04_resequenced_combined)){ 
    if(grepl(as.character(duplicated_file_data$read_id_list[i]), as.character(r04_resequenced_combined$read_id[j]))){ 
      #print(as.character(duplicated_file_data$read_id_list[i]))
      #print(as.character(r04_resequenced_combined$read_id[j]))
      duplicates_to_keep<-append(duplicates_to_keep, j)
        }
      }
    }
  }

r04_resequenced_combined<-r04_resequenced_combined[duplicates_to_keep,]



r04_isizes<-r04_resequenced_combined$isize
r04_sample<-c(rep('121r04_resequenced',length(r04_isizes)))
temp_table<-table(r04_isizes)
r04_freq_df<-data.frame(temp_table)
colnames(r04_freq_df)[1]<-'isize'
colnames(r04_freq_df)[2]<-'freq'
r04_freq_df$sample<-('121r04_resequenced')
r04_freq_df$percent<-100*r04_freq_df$freq/sum(r04_freq_df$freq)

human_isizes<-median_list
human_sample<-c(rep('Human',length(median_list)))
temp_table<-table(median_list)
human_freq_df<-data.frame(temp_table)
#human_freq_df<-data.frame(human_isizes,human_sample)
colnames(human_freq_df)[1]<-'isize'
colnames(human_freq_df)[2]<-'freq'
human_freq_df$sample<-('Human')
human_freq_df$percent<-100*human_freq_df$freq/sum(human_freq_df$freq)

freq_df<-rbind(human_freq_df,r04_freq_df)
freq_df$isize<-as.numeric(freq_df$isize)
freq_df$freq<-as.numeric(freq_df$freq)

head_freq_df<-freq_df[1:100000,]

####
#human_freq_df<-freq_df[freq_df$sample=='Human',]
percent_plot<-ggplot(human_freq_df, aes(x= human_freq_df$isize, y=human_freq_df$percent, color = human_freq_df$sample)) + 
  geom_point() + 
  xlab('Fragment length') + 
  ylab('percent of mapped reads in the sample') + 
  #xlim(c(0,500)) + 
  geom_vline(xintercept = 130) + 
  theme_classic()
percent_plot
ggsave(plot = percent_plot,'human_small.pdf')


test<-c()
d2 <- rep(human_freq_df$isize, human_freq_df$freq)
for(i in 1:nrow(human_freq_df)){ 
  d2 <- rep(d$Score, d$Frequency)

    }


####


percent_plot<-ggplot(head_freq_df, aes(x= head_freq_df$isize, y=head_freq_df$percent, color = head_freq_df$sample)) + 
  geom_line() + 
  xlab('Fragment length') + 
  ylab('percent of mapped reads in the sample') + 
  xlim(c(50,400)) + 
  geom_vline(xintercept = 167) +
  theme_classic()
percent_plot



percent_plot<-ggplot(freq_df, aes(x= freq_df$isize, y=freq_df$percent, color = freq_df$sample)) + 
  geom_line() + 
  xlab('Fragment length') + 
  ylab('% of mapped reads in the sample') +
  xlim(c(0,500))+
  geom_vline(xintercept = 168) + 
  theme_classic()
percent_plot

ggsave('121_r04_only_realigned_human_vs_resequenced_frequency_plot.pdf', percent_plot)


#data for percent plot
write.csv(freq_df,'121_r04_only_realigned_human_vs_resequenced_frequency_plot_data.csv')




###calculating true human mean fragment length from resequenced sample 
resequenced_isizes<-read.csv('/Volumes/Seagate8Tb1/resquenced_121R04_D01_CFFv1_NB0222/human_fragment_lengths.txt')
mean_resequenced_human_fragment_length<-mean(resequenced_isizes$X0)
median_resequenced_human_fragment_length<-median(resequenced_isizes$X0)

mean_resequenced_cmv_fragment_length<-mean(r04_resequenced_combined$isize)
median_resequenced_cmv_fragment_length<-median(r04_resequenced_combined$isize)


pvalues<-c()
for(i in 1:50000){ 
  progress(i,50000)
  temp_cmv<-sample(r04_resequenced_combined$isize,3000)
  temp_human<-sample(resequenced_isizes$X0,3000)
  #result<-cor.test(temp_cmv,temp_human, method = "pearson")
  result<-wilcox.test(temp_human,temp_cmv)
  pvalues<-append(result$p.value,pvalues)
  }
plot(pvalues)



#p value on insert sizes 
t.test(human_isizes,r04_isizes, alternative = "two.sided")

tests<-c()
for(i in 1:100000){ 
  progress(i,100000)
  temp<-(sample(human_isizes, length(r04_isizes)))
  #print((temp))
  ttest<-(t.test(temp,r04_isizes, paired = TRUE))
  tests<-append(tests,ttest$p.value)
  }

testsdf<-data.frame(tests)

ttest_graph<-ggplot(testsdf, aes(x = testsdf$tests)) + 
  geom_freqpoly(bins = 1000) +
  xlim(c(0,summary(testsdf$tests)[5]))+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  ylim(0,2000) + 
  xlab('p-value') + 
  ylab('frequency (100000 iterations)') + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
 # ylim(c(0,8))
ttest_graph
ggsave(plot =ttest_graph, 'insert_size_p-value_distribution.pdf')

#library('dyplr')
#boxplot(tests)


human_percentile<-ecdf(human_isizes)
cmv_med_percentile<-human_percentile(median(r04_isizes))


tests<-c()
for(i in 1:10000){ 
  progress(i,10000)
  temp<-(sample(human_isizes, length(r04_isizes)))
  #print((temp))
  fktest<-(fligner.test(temp,r04_isizes, paired = TRUE))
  tests<-append(tests,fktest$p.value)
}

fk_df<-data.frame(tests)

fktest_graph<-ggplot(fk_df, aes(x = fk_df$tests)) + 
  geom_freqpoly(bins = 1000) +
  xlim(c(0,summary(fk_df$tests)[5]))+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  ylim(0,2000) + 
  xlab('p-value') + 
  ylab('frequency (100000 iterations)') + 
  theme_classic() + 
  theme(axisw.text.x = element_text(angle = 90, hjust = 1)) 
# ylim(c(0,8))
fktest_graph 
ggsave(plot =fktest_graph, 'insert_size_variance_distribution.pdf')






