library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(svMisc)
library(xlsx)
library(Rsamtools)


setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download')

qpcr_results<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/qpcr_results.xls', sheetIndex =3)
sample_read_data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/all_deduplicated/all_sample_data.csv')

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


data<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/qpcr_graph_update/6131_cmv_percent_vs_qpcr_load.xlsx', sheetIndex = 1)

#deduplicated read counts version
data<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/qpcr_graph_update/82419_deduplicated_cmv_percent_vs_qpcr_load.xlsx', sheetIndex = 1)


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


fit1 <- lm((data$CMv_quantity.ml_plasma) ~ (data$rpm ), data = data)
library(ggpmisc)

#with rpm
cmv_percent_vs_viral_laod_regression<-ggplot(data, aes(x = data$CMv_quantity.ml_plasma, y= data$rpm)) + 
  geom_point(aes( color = data$classification.1)) + 
  theme_classic() + 
  theme(legend.position="none") + 
  scale_y_log10() + 
  scale_x_log10(limits = c(NA, 1000)) +
  xlab('CMV copies/ml by qPCR') + 
  ylab('RPM of total reads mapping to CMV')  + 
  geom_smooth(method = "lm", se = FALSE)   + 
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), formula = fit1, parse = TRUE, size = 5) 

  #geom_abline(slope = fit1$coefficients[2], intercept = fit1$coefficients[1])
cmv_percent_vs_viral_laod_regression
setwd('/Users/gerbix/Documents/vikas/NIPT/all_deduplicated')
ggsave('rpm reads mapping to cmv with regression line_v2.pdf' ,cmv_percent_vs_viral_laod_regression, width = 3, height = 3)




