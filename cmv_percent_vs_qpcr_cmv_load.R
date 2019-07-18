library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)
library(xlsx)
library(Rsamtools)


setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches')


#qpcr_results<-read.xlsx2('/Users/gerbix/Documents/vikas/NIPT/31119_download/qpcr_results.xls', sheetIndex =3, colClasses=rep("character",9))
qpcr_results <- read_excel("qpcr_results.xls", 
                           sheet = "combined", col_types = c("text", 
                                                             "text", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "text", "text"))
sample_read_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/all_sample_data.csv')

combined<-qpcr_results[qpcr_results$`CMV_Quantity(10UL DNA)`> 0,]
combined_names<-as.character(combined$`Sample Name`)


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


data<-read.xlsx2('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/cmv_percent_vs_qpcr_load.xlsx', sheetIndex = 1)



data<-read_excel("/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/6131_cmv_percent_vs_qpcr_load.xlsx", 
                      col_types = c("text", "numeric", "numeric", 
                                              "numeric", "numeric", "numeric", 
                                              "numeric", "text", "numeric", "text", 
                                              "numeric", "numeric", "numeric", 
                                              "numeric", "text"))


data$percent<-data$count/data$read_counts

cmv_percent_vs_viral_laod<-ggplot(data, aes(x = data$`CMv_quantity/ml_plasma`, y= data$percent, color = data$classification...15)) + 
  geom_point() + 
  theme_classic() + 
  theme(legend.position="none") + 
  scale_y_log10() + 
  scale_x_log10() + 
  xlab('CMV copies/ml by qPCR') + 
  ylab('Percent of total reads mapping to CMV')  
cmv_percent_vs_viral_laod

ggsave('61319_cmv_percent_vs_viral_load.pdf', cmv_percent_vs_viral_laod, width = 3, height = 3)


#testing
cmv_percent_vs_viral_laod_regression<-ggplot(data, aes(x = data$`CMv_quantity/ml_plasma`, y= data$percent)) + 
  geom_point(aes( color = data$classification...15)) + 
  theme_classic() + 
  theme(legend.position="none") + 
  scale_y_log10() + 
  scale_x_log10() + 
  xlab('CMV copies/ml by qPCR') + 
  ylab('Percent of total reads mapping to CMV')  + 
  geom_smooth(method = "lm", se = FALSE)
cmv_percent_vs_viral_laod_regression
ggsave('61319_percent of reads mapping to cmv with regression line.pdf' ,cmv_percent_vs_viral_laod_regression, width = 3.5, height = 3.5)

cmv_percent_vs_viral_laod_regression<-ggplot(data, aes(x = data$`CMv_quantity/ml_plasma`, y= data$rpm)) + 
  geom_point(aes( color = data$classification...15)) + 
  theme_classic() + 
  theme(legend.position="none") + 
  scale_y_log10() + 
  scale_x_log10() + 
  xlab('CMV copies/ml by qPCR') + 
  ylab('RPM of total reads mapping to CMV')  + 
  geom_smooth(method = "lm", se = FALSE)
cmv_percent_vs_viral_laod_regression
ggsave('61319_RPM of reads mapping to cmv with regression line.pdf' ,cmv_percent_vs_viral_laod_regression, width = 3.5, height = 3.5)

fit1 <- lm(data$percent ~  data$`CMv_quantity/ml_plasma`)







#coverage graph: 
verified_reads<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/blast_positive_sequence_info.csv')
ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

verified_reads<-ReadFasta('/Users/gerbix/Documents/vikas/NIPT/31119_download/resequenced/resequenced_repeatmasked_7_removed.fasta')
bam_files<-list.files('/Users/gerbix/Documents/vikas/NIPT/31119_download/cmv_download/bams', pattern = '*.bam$')

verified_reads$name<-as.character(verified_reads$name)
#verified_reads$read_id<-NA
#verified_reads$sample_id<-NA

verified_reads$base_read_id <- NA
for(i in 1:nrow(verified_reads)){
  verified_reads$base_read_id[i]<-strsplit(verified_reads$name[i],'/')[[1]][1]
  #verified_reads$sample_id[i]<-strsplit(verified_reads$nameslist[i],'[.]')[[1]][1]
  }

verified_reads$pos<-NA
unique_samples<-unique(verified_reads$base_read_id)
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
  geom_freqpoly(bins = 30) + 
  scale_y_log10() + 
  #ylim() +  
  xlab('position') + 
  theme_classic() 
coverage_plot

ggsave('30_bincmv_coverage_plot.pdf', coverage_plot,width = 3, height = 3)





all_pos<-c()
for(i in 1:nrow(verified_reads)){ 
  progress(i, nrow(verified_reads))
  for(j in 0:36) { 
    temp<-(verified_reads$pos[i] + j)
    all_pos<-append(all_pos, temp)
    
    }
  }

positions<-data.frame(all_pos)
positions_nozero<-all_pos[all_pos > 0]
positions_over_one<-data.frame(positions_nozero)


coverage_plot<-ggplot(positions_over_one, aes ( x = positions_over_one$positions_nozero)) + 
  geom_freqpoly(bins= 50) + 
  scale_y_log10() + 
  theme_classic()
coverage_plot








