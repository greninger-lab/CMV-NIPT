#load fig 3 data
setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/maternal')

cmv<-combined[combined$type=="CMV",]
cmv$isize<-as.numeric(cmv$isize)

write.csv(cmv,'121r04_cmv_cumulative_dist.csv')

#load 4b data
setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/CMV/SOT')

cmv_sot<-combined_isize_df[combined_isize_df$type=='CMV',]
write.csv(cmv_sot,'SOT_cmv_cumulative_dist.csv')
