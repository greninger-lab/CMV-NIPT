setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal/fragment_patch')

strong_positive_sample_names<-df_for_graph$sample[df_for_graph$classification == 'strong positive']

bam_folder<-'/Users/gerbix/Documents/vikas/NIPT/all_deduplicated/original_bams_deduplicated/'

blasted_reads<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal/blast_positive_sequence_info.csv')

all_sample_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal/all_sample_data.csv')

blasted_reads$sample<-NA
blasted_reads$unique_id<-NA
for(i in 1:nrow(blasted_reads)){ 
  blasted_reads$sample[i]<-strsplit(as.character(blasted_reads$nameslist[i]),'[.]')[[1]][1]
  blasted_reads$unique_id[i]<-strsplit(as.character(blasted_reads$nameslist[i]),'-')[[1]][2]
}

to_keep<-as.character(all_sample_data$sample[all_sample_data$classification=='strong positive'])
blasted_reads<-blasted_reads[which(blasted_reads$sample %in% to_keep),]


blasted_reads$isize<-NA
for(i in 1:nrow(blasted_reads)){ 
  print(i)
  temp_dir<-paste0(bam_folder,'/',blasted_reads$sample[i], '.sam.bam.sorted.bam_dedup.bam')
  temp_bam<-scanBam(temp_dir)
  blasted_reads$isize[i]<-abs(temp_bam[[1]]$isize[which(grepl(blasted_reads$unique_id[i],temp_bam[[1]]$qname))])
}
blasted_reads<-blasted_reads[complete.cases(blasted_reads$isize),]
r04_removed<-blasted_reads[-which(blasted_reads$sample == '121R04_D01_CFFv1_NA0144'),]

summary(r04_removed$isize[blasted_reads$isize < 500], na.omit= TRUE)
