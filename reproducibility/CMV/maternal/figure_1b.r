library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(xlsx)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal')

int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] 

#Modify for your own all_sample_data.csv 
df_for_graph<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal/fragment_patch/all_sample_data_og.csv')
figure_1b_plot_fpm<-ggplot(df_for_graph, aes(x=fpm, fill = classification)) + 
   geom_histogram(bins = 45, color = 'black') +  
  scale_x_log10() + 
  theme_classic() +
  scale_fill_manual(labels = c( "FPM < 0.3", "FPM > 0.3"),values = c('#BEBADA',"#E27184")) +
  theme(legend.title=element_blank())+ 
  theme(legend.position = c(0.8, 0.8))+
  xlab('FPM') + 
  ylab('Number of samples') + 
  theme(text = element_text(size=10))+ 
  scale_y_continuous(expand = c(0,0), breaks = int_breaks)
figure_1b_plot_fpm
ggsave(plot=figure_1b_plot_fpm, filename = 'figure_1b.pdf', height= 4, width = 4, useDingbats=FALSE)

