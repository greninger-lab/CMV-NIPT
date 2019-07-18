library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(xlsx)

setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches')
qpcr_results<-read.xlsx2('qpcr_results.xls', sheetIndex = 3)
#read_counts <- read.csv("/Users/gerbix/Documents/vikas/NIPT/21419_download/bam_fragment_read_counts.csv")
cmv_all_positives<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/all_sample_data.csv', colClasses=c('character', 'character', 'character','character', 'character', 'character','character')) 



qpcrplot_copies_microliter<-ggplot(qpcr_results, aes(x=qpcr_results$BGLOBIN_Quantity.10UL.DNA., y=qpcr_results$CMV_Quantity.10UL.DNA., color=qpcr_results$classification)) + 
  theme_classic() +
  theme(legend.position= "bottom")+
  xlab('Beta-globin copies per reaction') +
  ylab('CMV copies per reaction') +
  ylim(c(0,100))+
  xlim(c(0,2500))+
#  scale_y_continuous(expand = c(0,0)) + 
#  scale_x_continuous(expand = c(0,0)) + 
  labs(color = 'Classification by RPM and number of reads') +
  geom_point(size = qpcr_results$size)
qpcrplot_copies_microliter

ggsave(plot = qpcrplot_copies_microliter, filename = 'qpcr_graphs/qpcr_cmv_bglobin_quants_per_10ul.pdf', height = 7, width = 7, useDingbats=FALSE)


qpcrplot_copies_ml_plasma<-ggplot(qpcr_results, aes(x=qpcr_results$bglobin_quantity.ml_plasma, y=qpcr_results$CMv_quantity.ml_plasma, color=qpcr_results$classification)) + 
  theme_classic() +
  theme(legend.position= "bottom")+
  xlab('Beta-globin copies per mL Plasma') +
  ylab('CMV copies per mL plasma') +
  #ylim(c(0,100))+
  #xlim(c(0,100))+
    scale_y_log10() + 
    scale_x_log10() + 
  labs(color = 'Classification by RPM and number of reads') +
  geom_point(size = 3 )
qpcrplot_copies_ml_plasma

ggsave(plot = qpcrplot_copies_ml_plasma, filename = 'qpcr_cmv_bglobin_quants_per_ml_plasma.pdf', height = 7, width = 7, useDingbats=FALSE)

read_counts <- read.csv("~/Documents/vikas/NIPT/clip_removed/cmv_full/read_counts.csv")


read_counts$cmv_bglobin_ratio_ml<-NA
for(i in 1:nrow(qpcr_results)){ 
  for(j in 1:nrow(read_counts)){
  if(grepl(as.character(qpcr_results$Sample.Name[i]), as.character(read_counts$sample[j]))){ 
    read_counts$cmv_bglobin_ratio_ml[j]<-qpcr_results$CMv_quantity.ml_plasma[i]/qpcr_results$bglobin_quantity.ml_plasma[i]
    
      }
    }
  }





ggplot(cmv_all_positives, aes(x=cmv_all_positives$read_counts) + 
  geom_histogram() + 
  scale_x_log10(breaks=c(1,2,3,4,5,6,10,20,50,100,200,500,1000))) + 
  geom_vline(aes(xintercept=15),color="blue", linetype="dashed", size=1) + 
  theme(panel.background = element_blank())






#figure 1b 

#controls<-(read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/controls_list.csv'))
#to_remove<-which(controls$x %in% df_for_graph$sample)


figure_1b_plot<-ggplot(cmv_all_positives, aes(x=read_counts)) + 
  geom_histogram(bins = 20, color= 'black') + 
  scale_x_log10() + 
  theme_classic() +
  xlab('number of reads') + 
  ylab('Number of samples containing this many reads')+ 
  scale_y_continuous(expand = c(0,0))
  #geom_vline(aes(xintercept=15),color="blue", linetype="dashed", size=1) 
figure_1b_plot
ggsave(plot = figure_1b_plot, filename = 'qpcr_graphs/figure_1b_plot.pdf', height = 4, width = 4, useDingbats = FALSE)


int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] 


df_for_graph<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/34_mismatches/all_sample_data.csv')
figure_1b_plot_rpm<-ggplot(df_for_graph, aes(x=rpm, fill = classification)) + 
   geom_histogram(bins = 45, color = 'black') +  
  scale_x_log10() + 
  #stat_bin(geom="text", colour="white", size=3.5, aes(label=)) +
  theme_classic() +
  theme(legend.position='none') + 
  xlab('RPM of Fragments') + 
  ylab('Number of samples containing this RPM') + 
  scale_y_continuous(expand = c(0,0), breaks = int_breaks)
figure_1b_plot_rpm
ggsave(plot=figure_1b_plot_rpm, filename = 'figure_1b_fragment_rpm_plot.pdf', height= 4, width = 4, useDingbats=FALSE)




figure_1b_plot_rpm<-ggplot(df_for_graph, aes(x=count, fill = classification)) + 
  geom_histogram(bins = 45, color='black') +  
  scale_x_log10() + 
  #scale_y_log10() +
  #stat_bin(geom="text", colour="white", size=3.5, aes(label=)) +
  theme_classic() +
  #theme(legend.position='none') + 
  xlab('Number of fragments') + 
  ylab('Number of samples containing this RPM') 
  #scale_y_continuous(expand = c(0,0))
figure_1b_plot_rpm
ggsave(plot=figure_1b_plot_rpm, filename = 'qpcr_graphs/number_of_fragments.pdf', height= 4, width = 4, useDingbats=FALSE)








