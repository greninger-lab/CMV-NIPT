library(ggplot2)


resequenced_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/resequenced_human_insert_sizes.csv')
resequenced_data$run<-'resequenced'
original_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/original_human_insert_sizes.csv')
original_data$run<-'original'

all_reads<-rbind(resequenced_data[resequenced_data$type=='with_duplicates',], original_data[original_data$type=='with_duplicates',])
deduplicated_reads<-rbind(resequenced_data[resequenced_data$type=='deduplicated',], original_data[original_data$type=='deduplicated',])
duplicate_reads<-rbind(resequenced_data[resequenced_data$type=='extracted_duplicates',], original_data[original_data$type=='extracted_duplicates',])

#all of the reads
all_reads_plot<-ggplot(all_reads, aes( x = all_reads$isize , y = all_reads$percent, color = all_reads$run)) + 
  geom_vline(xintercept = 168) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line(aes(alpha=0.4)) 
all_reads_plot
ggsave(plot = all_reads_plot, 'all_reads_overlayed.pdf', height = 8, width = 8)



#deduplicated reads
deduplicated_reads_plot<-ggplot(all_reads, aes( x = deduplicated_reads$isize , y = deduplicated_reads$percent, color = deduplicated_reads$run)) + 
  geom_vline(xintercept = 168) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line(aes(alpha=0.4)) 
deduplicated_reads_plot
ggsave(plot = deduplicated_reads_plot, 'deduplicated_overlayed.pdf', height = 8, width = 8)

#duplicate reads
duplicate_reads_plot<-ggplot(duplicate_reads, aes( x = duplicate_reads$isize , y = duplicate_reads$percent, color = duplicate_reads$run)) + 
  geom_vline(xintercept = 168) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line(aes(alpha=0.4)) 
duplicate_reads_plot
ggsave(plot = duplicate_reads_plot, 'duplicates_overlayed.pdf', height = 8, width = 8)




duplicate_vs_non_duplicate_resequenced<-rbind(resequenced_data[resequenced_data$run=='resequenced' & resequenced_data$type=='extracted_duplicates',], resequenced_data[resequenced_data$run=='resequenced' & resequenced_data$type=='deduplicated',])

duplicate_vs_non_duplicate_resequenced_plot<-ggplot(duplicate_vs_non_duplicate_resequenced, aes( x = duplicate_vs_non_duplicate_resequenced$isize , y = duplicate_vs_non_duplicate_resequenced$percent, color = duplicate_vs_non_duplicate_resequenced$type)) + 
  geom_vline(xintercept = 168) + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line(aes(alpha=0.4)) 
duplicate_vs_non_duplicate_resequenced_plot






