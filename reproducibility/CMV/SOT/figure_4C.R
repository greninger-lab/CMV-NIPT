library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)
library(wesanderson)
library(RColorBrewer)
library(Hmisc)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT')
blast_hits_file<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/blast_hits.csv', header = FALSE, col.names = c('full_ID','count'))
filenames = list.files('/Users/gerbix/Documents/vikas/NIPT/new_samples',pattern = '.bam$', full.names = TRUE)
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
  temp_bam<-scanBam(filenames[i])
  if((identical(temp_bam[[1]]$qname, character(0)))){ 
    print(i)
    next
  }
  all_bam_read_IDs<-temp_bam[[1]]$qname  
  matches<-c()
  for(j in 1:nrow(blast_hits_file)){ 
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
  ylab('occurences')+
  xlab('fragment length')+
  annotate("text", x = 400, y = 15, label =  paste0('median=', median(combined$isize)))  +
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")
cmv_plot

combined_trimmed<-combined[combined$isize<501,]

cmv_cumulative_frequency<-ggplot(combined_trimmed, aes(x = combined_trimmed$isize, color = combined_trimmed$sample)) + 
  stat_ecdf(geom = 'step', size  =.5 ) +
  theme_classic() +  
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.position='bottom') + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') 
cmv_cumulative_frequency

#pulling human insert sizes from sample 244P16_H02_CFFv2_NB0289
human<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/p16.results.txt', sep = '', header = FALSE, col.names = c('frequency', 'isize'))
human<-human[human$isize>1,]

human_isizes<-rep(human$isize, human$frequency)
human_df<-as.data.frame(human_isizes)
colnames(human_df)[1]<-'isize'
human_df$sample<-'P16_Human'
human_df$read_id<-'eh'


human_df<-human_df[human_df$isize < 1000,]
human_cmv_combined<-rbind(human_df,combined_trimmed)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(8, "Set3"))
colourCount = length(unique(human_cmv_combined$sample))


cmv_cumulative_frequency<-ggplot(combined_trimmed, aes(x = isize, color = sample)) + 
  stat_ecdf(geom = 'step', size  =.5 )+
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_classic() +  
  theme(legend.position='bottom') + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') 
cmv_cumulative_frequency


colors<-c("#8DD3C7", "#FFFFB3" ,"#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#05188B")
cumulative_freq_with_human<-ggplot(human_cmv_combined, aes(x = human_cmv_combined$isize, color = human_cmv_combined$sample)) + 
  theme_classic() +  
  theme(legend.position='none') + 
  scale_color_manual(values = colors) + 
  xlab('Insert size') + 
  ylab ('Cumulative frequency') + 
  stat_ecdf(geom = 'step', size  =.5 ) 
cumulative_freq_with_human
ggsave(plot = cumulative_freq_with_human, 'figure_4c.pdf', height = 3, width = 3)
save.image("~/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/figure_4c.rdata")





