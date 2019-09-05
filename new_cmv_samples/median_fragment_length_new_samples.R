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
  #print(100 * (i / length(filenames)))
  #print(i)
  temp_bam<-scanBam(filenames[i])
  if((identical(temp_bam[[1]]$qname, character(0)))){ 
    #print(i)
    next
  }
  all_bam_read_IDs<-temp_bam[[1]]$qname  
  matches<-c()
  for(j in 1:nrow(blast_hits_file)){ 
    print(100 * (j/nrow(blast_hits_file)))
    for(k in 1:length(all_bam_read_IDs)){
     # print(blast_hits_file$read_ID[j])
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
cmv_plot<-ggplot(combined, aes( x=combined$isize, fill = sample)) + 
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(combined$isize)) +
  xlim(0,500) + 
  #ylim(0,200) + 
  ylab('occurences')+
  xlab('fragment length')+
  #annotate("text", x = 400, y = 18 , label =  paste0('mean=', mean(combined$isize))) + 
  annotate("text", x = 400, y = 600, label =  paste0('median=', median(combined$isize[combined$isize<500])))  +
  labs(colour="file") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() + 
  theme(legend.position="right")
cmv_plot
ggsave(plot = cmv_plot, 'new_cmv_insert_size_colored.pdf', height = 5, width = 5)


#overlaying human 
human<-read.csv('/Volumes/Seagate8Tb1/244P_data/244P16_H02_CFFv2_NB0289.final.bam.results.txt', sep = ' ', header = FALSE, col.names = c('frequency', 'isize'))
human<-human[human$isize>0,]
human_isizes<-rep(human$isize, human$frequency)
human_df<-as.data.frame(human_isizes)
colnames(human_df)[1]<-'isize'
human_df$sample<-'P16_Human'
human_df$read_id<-'eh'


cmv_df<-as.data.frame(table(combined$isize[combined$isize<500]))
cmv_df$percent<-NA
for(i in 1:nrow(x)){ 
  cmv_df$percent[i] <- 100 * (cmv_df$Freq[i] / sum(cmv_df$Freq))
  }
cmv_df$type<-'CMV'

human_isize_df<-as.data.frame(table(human_df$isize[human_df$isize<500]))
human_isize_df$percent<-NA
for(i in 1:nrow(x)){ 
  human_isize_df$percent[i] <- 100 * (human_isize_df$Freq[i] / sum(human_isize_df$Freq))
}
human_isize_df$type<-'Human'

combined_isize_df<-rbind(human_isize_df, cmv_df)


cmv_plot<-ggplot(combined_isize_df, aes( x=as.numeric(as.character(combined_isize_df$Var1)), y = as.numeric(as.character(combined_isize_df$percent)),color = combined_isize_df$type)) +
  geom_line() + 
  #geom_vline(xintercept = median(combined$isize)) +
  xlim(0,500) + 
  ylim(0,3) + 
  ylab('percent')+
  xlab('insert size')+
  #annotate("text", x = 400, y = 18 , label =  paste0('mean=', mean(combined$isize))) + 
  #annotate("text", x = 400, y = 600, label =  paste0('median=', median(combined$isize[combined$isize<500])))  +
  #labs(colour="file") +
  #scale_y_continuous(expand = c(0,0)) +
  #theme(legend.title ="") + 
  theme_classic() +
  #theme(legend.title = element_blank()) 
  theme(legend.position = 'none') 
cmv_plot

ggsave(plot = cmv_plot, 'new_cmv_isize_distribution.pdf', height = 3, width = 3)

######
for(i in (unique(combined$sample))){ 
  print(i)
  print(median(combined$isize[combined$sample==i & combined$isize < 500]))
  print('')
    }


P17<-combined[combined$sample=='244P17_A03_CFFv2_NB0289.sam.bam' & combined$isize < 500,]
median(P17$isize[P17$isize<500])

#finding size of the two peaks 
median(combined$isize[combined$isize < 100])

median(combined$isize[combined$isize < 250 & combined$isize > 100])


