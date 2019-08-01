library(ggplot2)



repeatmasked_fasta<-read.fasta('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/resequenced_repeatmasked_7_removed.fasta')
fasta_IDs<-names(repeatmasked_fasta)

fasta_IDs_trimmed<-c()
for( i in 1:length(fasta_IDs)){
  progress(i,length(fasta_IDs))
  fasta_IDs_trimmed<-append(fasta_IDs_trimmed,strsplit(fasta_IDs, '/')[[i]][1])
} 
fasta_IDs_deduplicated<-fasta_IDs_trimmed[-(which(duplicated(fasta_IDs_trimmed)))]

#CMV_blast_hits_deduplicated<-CMV_blast_hits[-(which(duplicated(CMV_blast_hits$read_ID))),]


isize<-c()
read_id<-c()
sample<-c()

temp_bam<-scanBam('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/cmv_just_duplicates.bam')
all_bam_read_IDs<-temp_bam[[1]]$qname  
matches<-c()
for(j in 1:length(fasta_IDs_deduplicated)){ 
  print(100 * j / length(fasta_IDs_deduplicated))
  for(k in 1:length(all_bam_read_IDs)){
    if(grepl( fasta_IDs_deduplicated[j], all_bam_read_IDs[k])) { 
      isize<-append(isize, temp_bam[[1]]$isize[k])
      read_id<-append(read_id, paste(temp_bam[[1]]$qname[k],'.1'))
      sample<-append(sample, 'CMV')
    }
  }
}

CMV_matched<-data.frame(sample, read_id, isize)
CMV_matched<-CMV_matched[CMV_matched$isize>0,]
CMV_matched<-CMV_matched[CMV_matched$isize<500,]

CMV_matched<-CMV_matched[complete.cases(CMV_matched$isize),]

CMV_frequencies<-data.frame(table(CMV_matched$isize))
colnames(CMV_frequencies)<-(c('isize','freq'))
CMV_frequencies$percent<- 100 * ( as.numeric(as.character(CMV_frequencies$freq)) / (sum(CMV_frequencies$freq))) 
CMV_frequencies$type<-'CMV'

#calculating CMV median
CMV_isize_exanded<-rep(CMV_frequencies$isize, CMV_frequencies$freq) 
CMV_isize_exanded<-CMV_isize_exanded[as.numeric(as.character(CMV_isize_exanded)) < 500]
CMV_median<-median(as.numeric(as.character(CMV_isize_exanded)))
CMV_frequencies$type<-'duplicates_only'



write.csv(CMV_matched,('cmv_duplicates_read_data.csv'))

cmv_deduplicated_read_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/cmv_duplicates_removed_read_data.csv')
cmv_deduplicated<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/deduplicating/cmv_duplicates_removed_insert_sizes.csv')
cmv_deduplicated$type<-'deduplicated'
cmv_just_duplicates<-CMV_frequencies


cmv_deduplicated$X<-NULL
all_cmv<-rbind(cmv_deduplicated, CMV_frequencies)



CMV_plot<-ggplot(all_cmv, aes( x = as.numeric(all_cmv$isize) , y = all_cmv$percent, color = all_cmv$type)) + 
  geom_vline(xintercept = 143, color = 'blue') + 
  theme_classic() + 
  xlim(c(0,500)) +
  xlab('insert size') + 
  ylab('percent')  +
  geom_line() 
CMV_plot
ggsave(plot = CMV_plot, 'cmv_deduplicated_vs_duplicates.pdf', width = 8, height =8)

#median calculations
#deduplicated cmv
CMV_deduplicated_isize_exanded<-rep(cmv_deduplicated$isize, cmv_deduplicated$freq) 
CMV_deduplicated_isize_exanded<-CMV_deduplicated_isize_exanded[as.numeric(as.character(CMV_deduplicated_isize_exanded)) < 500]
CMV_deduplicated_median<-median(as.numeric(as.character(CMV_deduplicated_isize_exanded)))

#cmv duplicate reads 

CMV_duplicates_isize_exanded<-rep(CMV_frequencies$isize, CMV_frequencies$freq) 
CMV_duplicates_isize_exanded<-CMV_duplicates_isize_exanded[as.numeric(as.character(CMV_duplicates_isize_exanded)) < 500]
CMV_duplicates_median<-median(as.numeric(as.character(CMV_duplicates_isize_exanded)))






