library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(xlsx)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV')

int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] 

#Modify for your own all_sample_data.csv 
df_for_graph<-read.csv('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/all_sample_data.csv')
figure_1b_plot_fpm<-ggplot(df_for_graph, aes(x=fpm, fill = classification)) + 
   geom_histogram(bins = 45, color = 'black') +  
  scale_x_log10() + 
  theme_classic() +
  theme(legend.position='none') + 
  xlab('FPM of Fragments') + 
  ylab('Number of samples containing this FPM') + 
  scale_y_continuous(expand = c(0,0), breaks = int_breaks)
figure_1b_plot_fpm
ggsave(plot=figure_1b_plot_rpm, filename = 'figure_1b_fragment_fpm_plot.pdf', height= 4, width = 4, useDingbats=FALSE)

