<<<<<<< HEAD
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
  #scale_colour_brewer(palette = 'Set2') + 
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
=======
#matching container IDs to sequencing IDs and creating standard curves for cmv samples 
library("ggplot2")
library("xlsx")

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT')

#matching qpcr values to fpm values for SOT cmv samples
cmv_ids<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/new_cmv_pulled_ids.csv', sep = '\t', header = FALSE, col.names = c('sequencing_id', 'bid', 'lmid'))
cmv_ids$bid<-as.character(cmv_ids$bid)
for(i in 1:nrow(cmv_ids)){ 
  cmv_ids$bid[i]<-strsplit(as.character(cmv_ids$bid[i]), 'P')[[1]][2]
}

rpkm_values<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/deduplicated/all_sample_data.csv')
rpkm_values$sample<-as.character(rpkm_values$sample)
rpkm_values$rpm<-as.numeric(as.character(rpkm_values$rpm))
cmv_quants<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/new_samples/CMVPulled.xlsx', sheetIndex = 1)

rpkm_values$quant<-NA
rpkm_values$sample<-as.character(rpkm_values$sample)
for(i in 1:nrow(cmv_ids)){ 
  temp<-which(grepl(cmv_ids$sequencing_id[i],rpkm_values$sample))
  print(rpkm_values$sample[temp])
  temp_quant<-which(grepl(cmv_ids$bid[i], as.character(cmv_quants$Barcode.ID.)))
  rpkm_values$quant[temp]<-cmv_quants$result_num[temp_quant]
}

rpkm_values$quant_adjusted<-rpkm_values$quant * 4
curve<-ggplot(rpkm_values, aes(x = rpm, y = quant)) + 
  geom_point() +
  ylim(c(0,20000)) +
  geom_smooth(method = "lm", se = FALSE, alpha = .5) + 
  theme_classic() + 
  geom_vline(xintercept = .3, linetype = 'dotted')
curve


write.csv(rpkm_values, 'fpm_values_with_quants.csv')

#matching qpcr values to fpm values for NIPT samples
original_cmv<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/qpcr_graph_update/6131_cmv_percent_vs_qpcr_load.xlsx', sheetIndex = 1)
original_reformatted<- original_cmv[,c(1,4,11,12,13,14,15)]
colnames(original_reformatted)<-c('sample','quant','count','rpkm','read_counts','rpm','classification')
original_reformatted$time<-'original'
original_reformatted$quant_adjusted<-original_reformatted$quant

rpkm_values$X<-NULL
rpkm_values$time<-'new'

original_new_combined<-rbind(original_reformatted,rpkm_values)
original_new_combined$time<-as.character(original_new_combined$time)

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(8, "Set3"))
colourCount = length(unique(original_new_combined$time))

plot<-ggplot(original_new_combined, aes(x = rpm, y = quant_adjusted, color= time)) + 
  geom_point() +
  scale_y_log10(breaks = c( 1, 10, 100, 1000, 10000,100000))+ 
  scale_x_log10(limits = c(.01, 100)) + 
  ylab('CMV copies/mL') + 
  xlab('CMV FPM') + 
  geom_smooth(method = "lm", se = FALSE, alpha = .5, aes(group=1), color = 'black') + 
  scale_color_manual(values = c("#729AF2","#BF6FF7")) + 
  theme_classic()  + 
  theme(legend.title=element_blank(), legend.position = 'none') 
plot
ggsave(plot = plot, 'figure_4c.pdf', height = 3, width = 3)

summary(lm(rpm ~ quant_adjusted, data=original_new_combined))


>>>>>>> fc8ed0222b2c2448d148e1835a033131476adb6a


