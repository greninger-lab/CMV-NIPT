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
lb1 <- paste0("r^2= ", signif(r2$adj.r.squared, digits = 2))
plot<-ggplot(original_new_combined, aes(x = rpm, y = quant_adjusted, color= time)) + 
  geom_point() +
  scale_y_log10(breaks = c( 1, 10, 100, 1000, 10000,100000), labels = trans_format("log10", math_format(10^.x)))+ 
  scale_x_log10(limits = c(.01, 100), labels = trans_format("log10", math_format(10^.x)))+  
  ylab('CMV copies/mL') + 
  xlab('CMV FPM') + 
  geom_smooth(method = "lm", se = FALSE, alpha = .5, aes(group=1), color = 'black') + 
  scale_color_manual(values = c("#729AF2","#BF6FF7"), labels = c("Maternal", "Transplant")) + 
  theme_classic()  + 
  annotate("text", x = 25  , y =.5, label = lb1, parse = FALSE, vjust =1, size = 3 ) +
  theme(text = element_text(size=8),
        legend.title=element_blank(),
        legend.position = c(.8,.2),
        legend.spacing.y = unit(0.01, 'cm')) 
plot
ggsave(plot = plot, 'figure_4a_modified.pdf', height = 3, width = 3)

save.image("~/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT/fragment_patch/figure_4A.rdata")




