library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)
library(wesanderson)
library(RColorBrewer)
library(Hmisc)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch')
blast_hits_file<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/blast_hits.csv', header = FALSE, col.names = c('full_ID','count'))
filenames = list.files('/Users/gerbix/Documents/vikas/NIPT/new_samples',pattern = '.bam$', full.names = TRUE)
plotslist<-c()


blast_hits_file$full_ID<-as.character(blast_hits_file$full_ID)
for( i in 1:nrow(blast_hits_file)){
  blast_hits_file$read_ID[i]<-strsplit(blast_hits_file$full_ID, '-')[[i]][2]
} 

blast_hits_file<-blast_hits_file[blast_hits_file$count>75,]
to_remove<-c()
duplicated<-which(duplicated(blast_hits_file$read_ID))
for(i in 1:length(duplicated)){ 
  duplicates<-which(blast_hits_file$read_ID == blast_hits_file$read_ID[duplicated[i]])
  first = duplicates[1]
  second = duplicates[2]
  #print(length(duplicates))
  #print(length(duplicates))
  #print(i)
  print(duplicates)
  # if( length( duplicates == 2)){
  #   print(i)
  # }
  if(blast_hits_file$count[first] >= blast_hits_file$count[second] & !identical(duplicates, integer(0))){ 
    if(length(duplicates) > 2){ 
      to_remove<-append(to_remove, duplicates[1:length(duplicates)])
      next
    }
    if(!identical(duplicates, integer(0))){ 
      to_remove<-append(to_remove, second)
    }
    else{ 
      next}
  }
  else{ 
    to_remove<-append(to_remove,first)
    
  }
}

to_remove<-to_remove[complete.cases(to_remove)]
blast_hits_file<-blast_hits_file[-to_remove,]
blast_hits_file<-blast_hits_file[-c(which(nchar(blast_hits_file$read_ID) > 5)),]

isize<-c()
read_id<-c()
sample<-c()
for(i in 1:length(filenames)){ 
  temp_bam<-scanBam(filenames[i])
  if((identical(temp_bam[[1]]$qname, character(0)))){ 
    #print(i)
    next
  }
  all_bam_read_IDs<-temp_bam[[1]]$qname  
  matches<-c()
  for(j in 1:nrow(blast_hits_file)){ 
    for(k in 1:length(all_bam_read_IDs)){
      #print(blast_hits_file$read_ID[j])
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


human_df<-human_df[human_df$isize < 500,]
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


colors<-c("#8DD3C7", "#F4EB42" ,"#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#05188B")
cumulative_freq_with_human<-ggplot(human_cmv_combined, aes(x = human_cmv_combined$isize, color = human_cmv_combined$sample)) + 
  theme_classic() +  
  scale_color_manual(values = colors, labels= c('P12','P13','P14','P15','P16','P17','P18','P16 (Human)')) + 
  xlab('Fragment size') + 
  theme(legend.title=element_blank()) + 
  theme(text = element_text(size=8)) +
  theme(legend.text=element_text(size=6)) +
  xlim(c(0,500)) +
  #geom_vline(xintercept = 500) +
  theme(legend.position = c(0.8, 0.3)) + 
  ylab ('Cumulative frequency') + 
  theme(legend.key.size = unit(.3, "cm")) + 
  stat_ecdf(geom = 'step', size  =.5, pad = FALSE) 
cumulative_freq_with_human
ggsave(plot = cumulative_freq_with_human, 'figure_4c.pdf', height = 3, width = 3)
save.image("~/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch/figure_4c.rdata")




