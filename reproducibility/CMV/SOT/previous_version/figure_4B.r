library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)

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
    next
  }
  all_bam_read_IDs<-temp_bam[[1]]$qname  
  matches<-c()
  for(j in 1:nrow(blast_hits_file)){ 
    print(100 * (j/nrow(blast_hits_file)))
    for(k in 1:length(all_bam_read_IDs)){
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
  ylab('occurences')+
  xlab('fragment length')+
  annotate("text", x = 400, y = 600, label =  paste0('median=', median(combined$isize[combined$isize<500])))  +
  labs(colour="file") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() + 
  theme(legend.position="right")
cmv_plot


#overlaying human 
human<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/p16.results.txt', sep = '', header = FALSE, col.names = c('frequency', 'isize'))
human<-human[human$isize>0,]
human_isizes<-rep(human$isize, human$frequency)
human_df<-as.data.frame(human_isizes)
colnames(human_df)[1]<-'isize'
human_df$sample<-'P16_Human'
human_df$read_id<-'eh'


cmv_df<-as.data.frame(table(combined$isize[combined$isize<500]))
cmv_df$percent<-NA
for(i in 1:nrow(cmv_df)){ 
  cmv_df$percent[i] <- 100 * (cmv_df$Freq[i] / sum(cmv_df$Freq))
  }
cmv_df$type<-'CMV'

human_isize_df<-as.data.frame(table(human_df$isize[human_df$isize<500]))
human_isize_df$percent<-NA
for(i in 1:nrow(human_isize_df)){ 
  human_isize_df$percent[i] <- 100 * (human_isize_df$Freq[i] / sum(human_isize_df$Freq))
}
human_isize_df$type<-'Human'

combined_isize_df<-rbind(human_isize_df, cmv_df)

cmv_plot_4b<-ggplot(combined_isize_df, aes( x=as.numeric(as.character(combined_isize_df$Var1)), y = as.numeric(as.character(combined_isize_df$percent)),color = combined_isize_df$type)) +
  geom_line() + 
  xlim(0,500) + 
  ylim(0,2.5) + 
  scale_color_manual(values=c("#ED6464", "#05188B"), labels = c("CMV", "Human")) +
  theme(legend.position = c(0.8, 0.6)) + 
  ylab('percent')+
  xlab('Fragment size')+
  theme_classic() +
  theme(text = element_text(size=8)) +
  theme(legend.position = c(0.8, 0.6)) + 
  theme(legend.title=element_blank())
  cmv_plot_4b

ggsave(plot = cmv_plot_4b, 'figure_4B.pdf', height = 3, width = 3)
save.image("~/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/figure_4B.Rdata")
