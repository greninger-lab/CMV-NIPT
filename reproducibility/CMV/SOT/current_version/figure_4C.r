library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch')
blast_hits_file<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/deduplicated/blast_hits.csv', header = FALSE, col.names = c('full_ID','count'))
filenames = list.files('/Users/gerbix/Documents/vikas/NIPT/new_samples',pattern = '.bam$', full.names = TRUE)
plotslist<-c()


blast_hits_file$full_ID<-as.character(blast_hits_file$full_ID)
for( i in 1:nrow(blast_hits_file)){
  progress( i , nrow(blast_hits_file))
  blast_hits_file$read_ID[i]<-strsplit(blast_hits_file$full_ID, '-')[[i]][2]
} 

blast_hits_file$blast_pass <- TRUE
for(i in 1:nrow(blast_hits_file)){ 
  progress( i, nrow(blast_hits_file))
  if(blast_hits_file$count[i] > 4 ){ 
    blast_hits_file$blast_pass[i] <- TRUE
  }
  else{ 
    blast_hits_file$blast_pass[i]<-FALSE
  }
  
}

#blast_hits_file<-blast_hits_file[blast_hits_file$count>75,] 
for(i in 1:nrow(blast_hits_file)){ 
  if(blast_hits_file$blast_pass[i]==FALSE){ 
    temp_id<-blast_hits_file$read_ID[i]
    if(temp_id %in% blast_hits_file$read_ID[which(duplicated(blast_hits_file$read_ID))]){ 
      temp_indexes<-which(blast_hits_file$read_ID==temp_id)
      first<-temp_indexes[1]
      second<-temp_indexes[2]
      if(blast_hits_file$blast_pass[first] == TRUE | blast_hits_file$blast_pass[second]== TRUE ){ 
        blast_hits_file$blast_pass[first] <-TRUE
        blast_hits_file$blast_pass[second] <-TRUE
      }
    }
    
  }
}
#tf 13112:10034:9920 , 13610:3061:1660
#ff 21401:15559:5891 , 21509:19330:6011

blast_hits_file<-blast_hits_file[which(!(duplicated(blast_hits_file$read_ID))),]

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
ggsave(plot = cumulative_freq_with_human, 'figure_4C.pdf', height = 3, width = 3)
save.image("~/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch/figure_4c.rdata")




