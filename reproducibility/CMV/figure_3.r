library(ggplot2)
library(Rsamtools)
library(svMisc)
library(seqinr)
library(reshape2)
library(cowplot)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV')

# modify to include the resequenced 121r04 frequency table
filenames = '/Volumes/Seagate8Tb1/resquenced_121R04_D01_CFFv1_NB0222/aligned_to_hg38/duplicates_removed_read_lengths.txt'
plotslist<-c()


dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()

#extracts all insert sizes below 500
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
read_counts$sample<-as.character(read_counts$sample)

human_frequencies$percent<-(100 * human_frequencies$freq) / (sum(human_frequencies$freq))

#human median calculation
human_isize_expanded<-rep(human_frequencies$isize, human_frequencies$freq) 
human_isize_expanded_below250<-human_isize_expanded[human_isize_expanded < 500]
human_median<-median(as.numeric(as.character(human_isize_expanded_below250)))


human_frequencies$isize<-as.numeric(as.character(human_frequencies$isize))

#plot of just human insert sizes
human_plot<-ggplot(human_frequencies, aes( x = human_frequencies$isize , y = human_frequencies$percent)) + 
  geom_vline(xintercept = human_median) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
human_plot


#reading in CMV BAMs

#provide repeatmasked fasta for 121r04 resequenced
repeatmasked_fasta<-read.fasta('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/resequenced_repeatmasked_7_removed.fasta')
fasta_IDs<-names(repeatmasked_fasta)

fasta_IDs_trimmed<-c()
for( i in 1:length(fasta_IDs)){
  progress(i,length(fasta_IDs))
  fasta_IDs_trimmed<-append(fasta_IDs_trimmed,strsplit(fasta_IDs, '/')[[i]][1])
} 

#converts reads to fragments
fasta_IDs_deduplicated<-fasta_IDs_trimmed[-(which(duplicated(fasta_IDs_trimmed)))]




isize<-c()
read_id<-c()
sample<-c()

#reads in bam to find insert sizes matched by read ID from fasta
temp_bam<-scanBam('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/cmv_duplicates_removed.bam')
  all_bam_read_IDs<-temp_bam[[1]]$qname  
  matches<-c()
  for(j in 1:length(fasta_IDs_deduplicated)){ 
    print(100 * j / length(fasta_IDs_deduplicated))
    for(k in 1:length(all_bam_read_IDs)){
      if(grepl( fasta_IDs_deduplicated[j], all_bam_read_IDs[k])) { 
        isize<-append(isize, temp_bam[[1]]$isize[k])
        read_id<-append(read_id, paste(temp_bam[[1]]$qname[k],'.1'))
        sample<-append(sample, 'CMV')
      }
    }
  }

CMV_matched<-data.frame(sample, read_id, isize)
CMV_matched<-CMV_matched[CMV_matched$isize>0,]
CMV_matched<-CMV_matched[CMV_matched$isize<500,]

CMV_matched<-CMV_matched[complete.cases(CMV_matched$isize),]

CMV_frequencies<-data.frame(table(CMV_matched$isize))
colnames(CMV_frequencies)<-(c('isize','freq'))
CMV_frequencies$percent<- 100 * ( as.numeric(as.character(CMV_frequencies$freq)) / (sum(CMV_frequencies$freq))) 
CMV_frequencies$type<-'CMV'

#calculating CMV median
CMV_isize_exanded<-rep(CMV_frequencies$isize, CMV_frequencies$freq) 
CMV_isize_exanded<-CMV_isize_exanded[as.numeric(as.character(CMV_isize_exanded)) < 500]
CMV_median<-median(as.numeric(as.character(CMV_isize_exanded)))


# plot of just resequenced CMV insert sizes
CMV_plot<-ggplot(CMV_frequencies, aes( x = as.numeric(as.character(CMV_frequencies$isize)) , y = CMV_frequencies$percent)) + 
  geom_vline(xintercept = 143, color = 'blue') + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
CMV_plot

write.csv(CMV_matched, 'cmv_duplicates_removed_read_data.csv')
write.csv(CMV_frequencies,'cmv_duplicates_removed_insert_sizes.csv')

#combining human and CMV graphs 
combined<-rbind(human_frequencies,CMV_frequencies)

#red: ED6464
#blue: 05188B

combined_plot <- ggplot(combined, aes ( x = as.numeric(as.character(combined$isize)), y = as.numeric(as.character(combined$percent)), color = combined$type)) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  #theme_classic() + 
  xlim(c(0,500)) +
  ylim(c(0,2.5)) +
  scale_color_manual(values=c("#ED6464", "#05188B")) +
  #geom_vline(xintercept = 168) + 
  #geom_vline(xintercept = 125) + 
  xlab('Insert size') + 
  ylab('Percent within each alignment') + 
  geom_line(size = .75) 
combined_plot
ggsave(plot = combined_plot, 'deduplicated_figure_3a_500bp_1026.pdf', height = 3, width = 3 )






#p-value calculations 
CMV_isize_exanded<-as.numeric(as.character(CMV_isize_exanded))
human_isize_expanded<-as.numeric(as.character(human_isize_expanded))


#two sided t test on insert sizes 
t.test(human_isize_expanded,CMV_isize_exanded, alternative = "two.sided")

#running t tests 10000 iterations 
tests<-c()
for(i in 1:10000){ 
  progress(i,10000)
  temp<-(sample(human_isize_expanded, length(CMV_isize_exanded)))
  #print((temp))
  ttest<-(t.test(temp,CMV_isize_exanded, paired = TRUE))
  tests<-append(tests,ttest$p.value)
}

testsdf<-data.frame(tests)

#CMV is in 14th percentile of human insert size distribution
human_percentile<-ecdf(human_isize_expanded)
cmv_med_percentile<-human_percentile(median(CMV_isize_exanded)) * 100

#fligner-killeen test for homogeneity of variances between samples
tests<-c()
for(i in 1:10000){ 
  progress(i,10000)
  temp<-(sample(human_isize_expanded, length(CMV_isize_exanded)))
  fktest<-(fligner.test(temp,CMV_isize_exanded, paired = TRUE))
  tests<-append(tests,fktest$p.value)
}

fk_df<-data.frame(tests)

#Kolmogorov-Smirnov test

tests<-c()
for(i in 1:10000){ 
  progress(i,10000)
  temp<-(sample(human_isize_expanded, length(CMV_isize_exanded)))
  #print((temp))
  kstest<-(ks.test(temp,CMV_isize_exanded))
  tests<-append(tests,kstest$p.value)
}

ks_df<-data.frame(tests)

median(ks_df$tests)

kstest_graph<-ggplot(ks_df, aes(x = ks_df$tests)) + 
  geom_freqpoly(bins = 1000) +
  xlim(c(0,summary(ks_df$tests)[5]))+
  ylim(0,2000) + 
  xlab('p-value') + 
  ylab('frequency (100000 iterations)') + 
  theme_classic() 
kstest_graph 
ggsave(plot =fktest_graph, 'insert_size_variance_distribution.pdf')


plot(ecdf(CMV_isize_exanded))
plot(ecdf(human_isize_expanded))




CMV_cdf <- ecdf(CMV_isize_exanded)
human_cdf <- ecdf(human_isize_expanded)

CMV_cdf
ks<-ks.test(CMV_isize_exanded,human_isize_expanded)

human_isize_df<-data.frame(human_isize_expanded)
cmv_isize_df<-data.frame(CMV_isize_exanded)

plot(ecdf(CMV_isize_exanded)) 

plot(ecdf(human_isize_expanded),add = TRUE, col = 2 )


human_cdf_df<-data.frame()

human_subsampled<-data.frame(sample(human_isize_expanded, length(CMV_isize_exanded)))
colnames(human_subsampled)[1]<-'isizes'
human_subsampled$type = 'human'

cmv_isize_df$type = 'cmv'
colnames(cmv_isize_df)[1]<-'isizes'

subsampled_df<-rbind(human_subsampled, cmv_isize_df)



cum_frequency<-ggplot(subsampled_df, aes(x = subsampled_df$isizes, color = subsampled_df$type)) + 
  theme_classic() +  
  scale_color_manual(values=c("#ED6464", "#05188B")) +
  theme(legend.position='none') + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') + 
  stat_ecdf(geom = 'step', size  =.75 ) + 
  xlim(c(0,500))
cum_frequency
ggsave(plot = cum_frequency, 'cmv_deduplicated_cum_frequency_1026.pdf',width = 3, height = 3 )


shapiro.test(subsampled_df$isizes[subsampled_df$type=='human'])


#combined plot #cum_frequency
#red: ED6464
#blue: 05188B
plot_grid(combined_plot, cum_frequency, labels = c('A','B'))
ggsave(plot = last_plot(), height = 3, width = 6, filename = 'figure_3_grid.pdf')


