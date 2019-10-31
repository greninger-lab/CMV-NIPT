library(ggplot2)
library(Rsamtools)
library(svMisc)
library(seqinr)
library(reshape2)
library(cowplot)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6')

load(file='/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/figure_5A.rdata')
load(file='/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/figure_5B.rdata')
load(file='/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/figure_5C.rdata')

plot_grid(p5,plot, cum_frequency ,labels = c('A','B','C'), ncol = 3, align = 'v')
ggsave(plot = last_plot(), height = 3, width = 8, filename = 'figure_5.pdf')
