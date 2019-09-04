#matching container IDs to sequencing IDs and creating standard curves for cmv samples 
library("ggplot2")
library("readxl")


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
  