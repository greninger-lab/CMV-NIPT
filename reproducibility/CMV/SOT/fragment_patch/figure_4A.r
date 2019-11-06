#matching container IDs to sequencing IDs and creating standard curves for cmv samples 
library("ggplot2")
library("xlsx")

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch')

#matching qpcr values to fpm values for SOT cmv samples
cmv_ids<-read.csv('/Users/gerbix/Documents/vikas/NIPT/new_samples/new_cmv_pulled_ids.csv', sep = '\t', header = FALSE, col.names = c('sequencing_id', 'bid', 'lmid'))
cmv_ids$bid<-as.character(cmv_ids$bid)
for(i in 1:nrow(cmv_ids)){ 
  cmv_ids$bid[i]<-strsplit(as.character(cmv_ids$bid[i]), 'P')[[1]][2]
}

rpkm_values<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch/all_sample_data.csv')
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


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

r2=summary(lm(rpm ~ quant_adjusted, data=original_new_combined))
plot<-ggplot(original_new_combined, aes(x = rpm, y = quant_adjusted, color= time)) + 
  geom_point() +
  scale_y_log10(breaks = c( 1, 10, 100, 1000, 10000,100000), labels = scientific_10)+ 
  scale_x_log10(limits = c(.01, 100), labels =scientific_10)+  
  ylab('CMV copies/mL') + 
  xlab('CMV FPM') + 
  geom_smooth(method = "lm", se = FALSE, alpha = .5, aes(group=1), color = 'black') + 
  scale_color_manual(values = c("#729AF2","#BF6FF7"), labels = c("Maternal", "SOT")) + 
  theme_classic()  + 
  annotate("text", x = 30  , y =.20, label = paste0('R^2 = ',signif(r2$adj.r.squared,digits = 2))) +
  theme(text = element_text(size=8)) + 
  theme(legend.title=element_blank(), legend.position = c(.8,.2),legend.background=element_blank()) 
plot
ggsave(plot = plot, 'figure_4a.pdf', height = 3, width = 3)

save.image("~/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch/figure_4A.rdata")




