library(ggplot2)
library(Rsamtools)
library(svMisc)
library(seqinr)
library(reshape2)
library(cowplot)


setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal')

bam_folder<-'/Users/gerbix/Documents/vikas/NIPT/all_deduplicated/original_bams_deduplicated/'

blasted_reads<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal/blast_positive_sequence_info.csv')

all_sample_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal/all_sample_data.csv')

blasted_reads$sample<-NA
blasted_reads$unique_id<-NA
for(i in 1:nrow(blasted_reads)){ 
  blasted_reads$sample[i]<-strsplit(as.character(blasted_reads$nameslist[i]),'[.]')[[1]][1]
  blasted_reads$unique_id[i]<-strsplit(as.character(blasted_reads$nameslist[i]),'-')[[1]][2]
  }

to_remove<-as.character(all_sample_data$sample[all_sample_data$classification=='strong positive'])
blasted_reads<-blasted_reads[-which(blasted_reads$sample %in% to_remove),]


blasted_reads$isize<-NA
for(i in 1:nrow(blasted_reads)){ 
  print(i)
      temp_dir<-paste0(bam_folder,'/',blasted_reads$sample[i], '.sam.bam.sorted.bam_dedup.bam')
      temp_bam<-scanBam(temp_dir)
      blasted_reads$isize[i]<-abs(temp_bam[[1]]$isize[which(grepl(blasted_reads$unique_id[i],temp_bam[[1]]$qname))])
}
blasted_reads<-blasted_reads[complete.cases(blasted_reads$isize),]

summary(blasted_reads$isize[blasted_reads$isize < 500])

length(unique(blasted_reads$unique_id))