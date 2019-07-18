library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)
library(xlsx)
library(Rsamtools)


setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download')

qpcr_results<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/qpcr_results.xls', sheetIndex =3)
sample_read_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/all_sample_data.csv')

combined<-qpcr_results[qpcr_results$CMV_Quantity.10UL.DNA.>0,]
combined_names<-as.character(combined$Sample.Name)


for(i in 1:length(combined_names)){ 
  if(grepl('-', combined_names[i])){ 
    combined_names[i]<-strsplit(combined_names[i], '-')[[1]][1]
    }
  }

for(i in 1:length(combined_names)) { 
  for(j in 1:length(sample_read_data)){ 
    if( grepl(combined_names[i], sample_read_data$sample[j])){ 
      print(combined_names[i])
      }
    }
}


data<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/cmv_percent_vs_qpcr_load.xlsx', sheetIndex = 1)


data$percent<-data$count/data$read_counts

cmv_percent_vs_viral_laod<-ggplot(data, aes(x = data$CMv_quantity.ml_plasma, y= data$percent, color = data$classification.1)) + 
  geom_point() + 
  theme_classic() + 
  theme(legend.position="none") + 
  scale_y_log10() + 
  scale_x_log10() + 
  xlab('CMV copies/ml by qPCR') + 
  ylab('Percent of total reads mapping to CMV')  
cmv_percent_vs_viral_laod

ggsave('cmv_percent_vs_viral_load.pdf', cmv_percent_vs_viral_laod, width = 3, height = 3)


#testing

cmv_percent_vs_viral_laod_regression<-ggplot(data, aes(x = data$CMv_quantity.ml_plasma, y= data$percent)) + 
  geom_point(aes( color = data$classification.1)) + 
  theme_classic() + 
  theme(legend.position="none") + 
  scale_y_log10() + 
  scale_x_log10() + 
  xlab('CMV copies/ml by qPCR') + 
  ylab('Percent of total reads mapping to CMV')  + 
  geom_smooth(method = "lm", se = FALSE)
cmv_percent_vs_viral_laod_regression
ggsave('percent reads mapping to cmv with regression line.pdf' ,cmv_percent_vs_viral_laod_regression, width = 3, height = 3)

fit1 <- lm(data$percent ~  data$CMv_quantity.ml_plasma)







#coverage graph: 
verified_reads<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/blast_positive_sequence_info.csv')

bam_files<-list.files('/Users/gerbix/Documents/vikas/NIPT/31119_download/cmv_download/bams', pattern = '*.bam$')

verified_reads$nameslist<-as.character(verified_reads$nameslist)
verified_reads$read_id<-NA
verified_reads$sample_id<-NA

for(i in 1:nrow(verified_reads)){ 
  verified_reads$read_id[i]<-strsplit(verified_reads$nameslist[i],'-')[[1]][2]
  verified_reads$sample_id[i]<-strsplit(verified_reads$nameslist[i],'[.]')[[1]][1]
  }

verified_reads$pos<-NA
unique_samples<-unique(verified_reads$sample_id)
for(i in 1:length(unique_samples)){ 
  progress(i, length(unique_samples))
  tempbamname<-paste0('/Users/gerbix/Documents/vikas/NIPT/31119_download/cmv_download/bams/',unique_samples[i],'.sam.bam')
  #print(tempbamname)
  tempbam<-scanBam(tempbamname)
  for(j in 1:nrow(verified_reads)){ 
    if(verified_reads$sample_id[j] == unique_samples[i]){ 
      read_index<-which(grepl(verified_reads$read_id[j],tempbam[[1]]$qname))
      verified_reads$pos[j]<-tempbam[[1]]$pos[read_index][1]
      #print(tempbam[[1]]$pos[read_index])
      }
    }
  }

coverage_plot<-ggplot(verified_reads, aes(x = verified_reads$pos)) + 
  geom_freqpoly(binwidth = 5) + 
  theme_classic() 
coverage_plot

ggsave('cmv_coverage_plot.pdf', coverage_plot,width = 3, height = 3)






