#matching container IDs to sequencing IDs and creating standard curves for cmv samples 
library("ggplot2")
library("xlsx")


cmv_ids<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/new_cmv_pulled_ids.csv', sep = '\t', header = FALSE, col.names = c('sequencing_id', 'bid', 'lmid'))
cmv_ids$bid<-as.character(cmv_ids$bid)
for(i in 1:nrow(cmv_ids)){ 
  cmv_ids$bid[i]<-strsplit(as.character(cmv_ids$bid[i]), 'P')[[1]][2]
}

rpkm_values<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/all_sample_data.csv')
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


curve<-ggplot(rpkm_values, aes(x = rpm, y = quant)) + 
  geom_point() +
#  scale_x_log10() + 
  ylim(c(0,20000)) +
  #xlim(c(0,3)) + 
  geom_smooth(method = "lm", se = FALSE, alpha = .5) + 
  theme_classic() + 
  geom_vline(xintercept = .3, linetype = 'dotted')
curve
ggsave( plot = curve, 'cmv_quants.pdf', height = 4, width = 4)

write.csv(rpkm_values, 'rpkm_values_with_quants.csv')
  






#Original CMV samples with qCPR data 

original_cmv<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/cmv_percent_vs_qpcr_load.xlsx', sheetIndex = 1)
original_reformatted<- original_cmv[,c(1,4,11,12,13,14,15)]
colnames(original_reformatted)<-c('sample','quant','count','rpkm','read_counts','rpm','classification')
original_reformatted$time<-'original'

rpkm_values$X<-NULL
rpkm_values$time<-'new'

original_new_combined<-rbind(original_reformatted,rpkm_values)
original_new_combined$time<-as.character(original_new_combined$time)

plot<-ggplot(original_new_combined, aes(x = rpm, y = quant, color= time)) + 
  geom_point() +
  scale_y_log10(limits=c(.1, 1000000)) + 
  scale_x_log10(limits = c(.01, 100)) + 
  geom_smooth(method = "lm", se = FALSE, alpha = .5, aes(group=1), color = 'black') + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.position = 'bottom') + 
  geom_vline(xintercept = .3, linetype = 'dotted')
plot
ggsave(plot = plot, 'cmv_original_new_quant_rpm.pdf', height = 3, width = 3)

summary(lm(rpm ~ quant, data=original_new_combined))




