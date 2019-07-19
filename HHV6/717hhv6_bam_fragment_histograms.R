#for icihhv6 in cfdna insert sizes

library(Rsamtools)
library(Biostrings)
library(seqinr)
library(ggplot2)

setwd('/Volumes/Seagate8Tb1/nipt/bams/bams_for_vikas/partial')
read_counts <- read.csv("~/Documents/vikas/NIPT/clip_removed/cmv_full/read_counts.csv")
filenames = list.files(pattern = '*.bam.results.txt$')
plotslist<-c()

savePlot <- function(myPlot, plotname) {
  pdf(plotname, width = 8.5, height = 11)
  print(myPlot)
  dev.off()
}

dflist<-c()
statslist<-c()
totallengths<-c()
totaloccurences<-c()
for (i in 1:length(filenames)){ 
  print(i)
  #print(filenames[i])
  if (file.size(filenames[i]) != 0) {
    histo<- read.table(filenames[i], quote="\"", comment.char="")
    colnames(histo)[2]<-'length'
    colnames(histo)[1]<-'occurences'
    histo$type<-strsplit(filenames[i],".sam")[[1]][1]
    histo<-histo[histo$length>0,]
    if(nrow(histo)>0){
      # for(j in 1:length(histo$occurences)){ 
      #   for(k in 1:histo$occurences[j]){ 
      #     #print(histo$length[i])
      #     statslist<-append(statslist, histo$length[j])      
      #   }
      # }
      # statslistremove<-which(statslist>900)
      # if(length(statslistremove)>0){ 
      #   statslist<-statslist[-statslistremove]
      # }
      toremove<-which(histo[,2]> 900)
      if(length(toremove)>0){ 
        histo<-histo[-toremove,]
      }
      totallengths= append(totallengths,histo$length)
      totaloccurences =append(totaloccurences,histo$occurences)
      # plot<-ggplot(histo, aes(y=histo$occurences, x=histo$length)) + 
      #   geom_point() + 
      #   geom_vline(xintercept = median(statslist)) +
      #   ylab('occurences')+
      #   xlab('length') +
      #   annotate("text", x = 400, y = max(histo$occurences)-5, label =  paste0('mean=', mean(statslist))) + 
      #   annotate("text", x = 400, y = max(histo$occurences)-6, label =  paste0('median=', median(statslist)))
      # 
      # print(plot)
      # #print(mean(statslist))
      # plotname = paste0(filenames[i],'.plot.pdf')
      # print(filenames[1])
      # print(i)
      # #ggsave(filename=plotname, device = pdf)
      # savePlot(plot, plotname)
      dflist[[i]]<-histo
    }
  }
}

########revised 3/04/19#################
combineddf <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()
read_counts$sample<-as.character(read_counts$sample)
#combineddf<-combineddf[combineddf$length>74,]
plot<-ggplot(combineddf, aes( x=combineddf$length, y = combineddf$occurences, color = combineddf$type)) + 
  geom_point() + 
  #geom_histogram(binwidth = 5) + 
  #  geom_histogram(binwidth = 5) +
  #geom_vline(xintercept = median(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom") 
  #annotate("text", x = 400, y = max(histo$occurences)-5, label =  paste0('mean=', mean(combineddf$length))) + 
  #annotate("text", x = 400, y = max(histo$occurences)-7, label =  paste0('median=', median(combineddf$length)))
plot
 ggsave('30219_revised_mean_median_human_bam_fragments.pdf', plot = last_plot())

#calculating mean human read length
temp_sum<-0
for(i in 1:nrow(combineddf)){ 
  temp_sum = temp_sum  + (combineddf$occurences[i] * combineddf$length[i])
  }
mean_human_read_length <- temp_sum / (sum(combineddf$occurences))


##############end revision##############

#calculating mean and median for human fragments
combineddf$totals<-NA
for(i in 1:nrow(combineddf)){ 
  combineddf$totals[i]<-combineddf$occurences[i]*combineddf$length[i]
  }
mean(combineddf$totals)

plot<-ggplot(combineddf, aes(y=combineddf$occurences, x=combineddf$length)) + 
  geom_point(aes(color = combineddf$type), size = .5) + 
  geom_vline(xintercept = mean(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  #labs(colour="file") + 
  theme_classic() + 
  theme(legend.position = 'none')
  #theme(legend.position="bottom")
plot
ggsave('occurences_all_bams_average_length.pdf',plot, ,height = 3, width = 3, useDingbats = FALSE)

savePlot(plot, 'all_bams_average_length.pdf')


combineddf$denominator_adjusted<-NA 


for (i in 1:length(read_counts$sample)){ 
  for( j in 1:length(combineddf$type)) { 
    temp<-strsplit(as.character(read_counts$sample[i]), split = 'final') [[1]][1]
    if( temp == combineddf$type[j]) { 
        combineddf$denominator_adjusted[j]= combineddf$occurences[j] / read_counts$count[i]      
      }
    }
  }




denominator_adjusted_plot<-ggplot(combineddf, aes(y=combineddf$denominator_adjusted, x=combineddf$length)) + 
  geom_point(aes(color = combineddf$type), size = 1) + 
  geom_vline(xintercept = mean(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurrences/total read count of sample')+
  xlab('length')+
  labs(colour="file") + 
  theme_classic() +
  theme(legend.position="none")

denominator_adjusted_plot
ggsave('denominator_adjusted_all_bams_average_length.pdf',denominator_adjusted_plot, ,height = 3, width = 3, useDingbats = FALSE)
savePlot(denominator_adjusted_plot, 'denominator_adjusted_all_bams_average_length.pdf')


#for rpkm
combineddf$rpkm<-NA
for (i in 1:length(read_counts$sample)){ 
  for( j in 1:length(combineddf$type)) { 
    temp<-strsplit(as.character(read_counts$sample[i]), split = 'final') [[1]][1]
    if( temp == combineddf$type[j]) { 
      #combineddf$rpkm[j]= (10e9*combineddf$occurences[j]) / (read_counts$count[i] * 3088286401) 
      combineddf$rpkm[j] =(combineddf$occurences[j])/((3088286401/1000)*(read_counts$count[i]/1e9))
    }
  }
}

rpkm_plot<-ggplot(combineddf, aes(y=combineddf$rpkm, x=combineddf$length)) + 
  geom_point(aes(color = combineddf$type), size = 1) + 
  geom_vline(xintercept = mean(combineddf$length)) +
  xlim(0,500) + 
  ylab('rpkm')+
  xlab('length')+
  labs(colour="file") + 
  theme_classic() +
  theme(legend.position="none")

rpkm_plot




#############################
#####HHV-6 positive files######
#############################
setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/host_filtered_hhv6/bams')
system('$PWD')
read_counts <- read.csv("/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/read_counts_all.csv")
system("for i in *.bam ; do echo processing $i ; samtools view -@ 8 $i | awk '{print $9}' | sort | uniq -c > $i.results.txt ; done")
#system('')
filenames = list.files(pattern = '*.bam.results.txt$')
plotslist<-c()

savePlot <- function(myPlot, plotname) {
  pdf(plotname, width = 8.5, height = 11)
  print(myPlot)
  dev.off()
}

totallengths<-c()
totaloccurences<-c()
dflist<-c()
toremove<-c()
statslist<-c()
statslistremove<-c()
for (i in 1:length(filenames)){ 
  #print(filenames[i])
  if (file.size(filenames[i]) != 0) {
  histo<- read.table(filenames[i], quote="\"", comment.char="")
  colnames(histo)[2]<-'length'
  colnames(histo)[1]<-'occurences'
  histo$type<-strsplit(filenames[i],".sam")[[1]][1]
  histo<-histo[histo$length>0,]
  if(nrow(histo)>0){
  for(j in 1:length(histo$occurences)){ 
    for(k in 1:histo$occurences[j]){ 
      #print(histo$length[i])
      statslist<-append(statslist, histo$length[j])      
    }
  }
  statslistremove<-which(statslist>900)
  if(length(statslistremove)>0){ 
    statslist<-statslist[-statslistremove]
  }
  toremove<-which(histo[,2]> 900)
  if(length(toremove)>0){ 
    histo<-histo[-toremove,]
  }
  totallengths= append(totallengths,histo$length)
  totaloccurences =append(totaloccurences,histo$occurences)
  plot<-ggplot(histo, aes(y=histo$occurences, x=histo$length)) + 
    geom_point() + 
    geom_vline(xintercept = median(statslist)) +
    ylab('occurences')+
    xlab('length') +
    annotate("text", x = 400, y = max(histo$occurences)-5, label =  paste0('mean=', mean(statslist))) + 
    annotate("text", x = 400, y = max(histo$occurences)-6, label =  paste0('median=', median(statslist)))
  
  print(plot)
  #print(mean(statslist))
  plotname = paste0(filenames[i],'.plot.pdf')
  print(filenames[1])
  print(i)
  #ggsave(filename=plotname, device = pdf)
  savePlot(plot, plotname)
  dflist[[i]]<-histo
    }
  }
}


# combined<-data.frame(totaloccurences, totallengths)
# 
# plot<-ggplot(combined, aes(y=combined$totaloccurences, x=combined$totallengths)) + 
#   geom_point() + 
#   geom_vline(xintercept = mean(combined$totallengths)) +
#   xlim(0,250) + 
#   ylab('occurences')+
#   xlab('length') 

####revised code as of 3/34/19######
combined <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()
for(i in 1:nrow(combined)){ 
  for(j in 1:combined$occurences[i]){ 
      lengthlist<-append(lengthlist,combined$length[i])
      samplelist<-append(samplelist,combined$type[i])
    }
}
combineddf<-data.frame(lengthlist,samplelist)

colnames(combineddf)[1]<-'length'
colnames(combineddf)[2]<-'sample'
#read_counts$sample<-as.character(read_counts$sample)
#combineddf<-combineddf[combineddf$length>74,]
plot<-ggplot(combineddf, aes( x=combineddf$length)) + 
  geom_histogram(binwidth = 5) + 
  #  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom") + 
  annotate("text", x = 400, y = 400, label =  paste0('mean=', mean(combineddf$length))) + 
  annotate("text", x = 400, y = 300, label =  paste0('median=', median(combineddf$length)))
plot
ggsave('30219_revised_mean_median_cmv_bam_fragments.pdf', plot = last_plot())
########end revision#########


read_counts$sample<-as.character(read_counts$sample)
for (i in 1:length(read_counts$sample)) {
  print(i)
  read_counts$sample[i]=print(strsplit(as.character(read_counts$sample[i]),split='.fastq.gz')[[1]][1])
}

read_counts$percentage<-NA
read_counts$summed<-NA
for(i in 1:length(read_counts$sample)){ 
    tempcount<-0
    for(j in 1:nrow(combineddf)){ 
      if(read_counts$sample[i]==combineddf$type[j]){ 
        print('ok')
        tempcount= tempcount + combineddf$occurences[j]
        }
    }
    read_counts$percentage[i] <-100*tempcount/(read_counts$count[i])
    read_counts$summed[i]<-tempcount
   # print(read_counts$percentage[i])
      }

#####
combineddf$denominator_adjusted<-NA 

for (i in 1:length(read_counts$sample)){ 
  for( j in 1:length(combineddf$type)) { 
    temp<-strsplit(as.character(read_counts$sample[i]), split = 'final') [[1]][1]
    if( temp == combineddf$type[j]) { 
      combineddf$denominator_adjusted[j]= combineddf$occurences[j] / read_counts$count[i]      
    }
  }
}

combineddf$percent<-0
for(i in 1:nrow(combineddf)){ 
  print(i)
  for(j in 1:nrow(frequencies)){ 
    if(combineddf$length[i] == frequencies$Var1[j]){ 
      combineddf$percent[i] <- frequencies$percent[j]
      }
    }
  }


plot<-ggplot(combineddf, aes( x=combineddf$length, y= combineddf$percent)) + 
  geom_histogram(binwidth = 5, aes(fill = combineddf$sample, y=(100 * ..count../sum(..count..)))) + 

  #geom_histogram(binwidth = 5, aes(fill=combineddf$sample)) + 
#  geom_histogram(binwidth = 5) +
  #geom_vline(xintercept = median(combineddf$length)) +
  xlim(0,500) + 
  ylab('percent')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")

plot
ggsave(plot = plot, 'hhv6_isizes_stacked.pdf', width = 8.5, height = 11)
#ggsave('all_cmv_positives_average_length.pdf',plot, height = 3, width = 3,useDingbats = FALSE)

write.csv(combined,'hhv_6_isizes.csv')





hhv6_low_cluster<-c('98P11_C02_CFFv1_NB0120','93R20_D03_CFFv1_NA0087','97P10_B02_CFFv1_NB0119','99P03_C01_CFFv1_NB0121','120R22_F03_CFFv1_NA0137','104P04_D01_CFFv1_NA0266')
hhv6_low_cluster<-as.character(hhv6_low_cluster)
combineddf$sample<-as.character(combineddf$sample)
combineddf$cluster<-'high cluster'
for(i in 1:nrow(combineddf)){ 
  for(j in 1:length(hhv6_low_cluster)){ 
    if(hhv6_low_cluster[j] %in% combineddf$sample[i]) { 
      combineddf$cluster[i]<-'low cluster'
      }
    }
  }

plot<-ggplot(combineddf, aes( x=combineddf$length)) + 
  geom_histogram(binwidth = 5, aes(fill=combineddf$cluster)) + 
  #  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")

plot


ggsave('filled_all_cmv_positives_average_length.pdf',plot, height = 3, width = 3,useDingbats = FALSE)

savePlot(plot, 'all_cmv_positives_average_length.pdf')



ggplot(combineddf, aes(y=combineddf$denominator_adjusted, x=combineddf$length)) + 
  geom_point() + 
  geom_vline(xintercept = mean(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurrences/total read count of sample')+
  xlab('length')+
  labs(colour="file") + 
  theme_classic() +
  theme(legend.position="none")



#for rpkm
combineddf$rpkm<-NA
for (i in 1:length(read_counts$sample)){ 
  for( j in 1:length(combineddf$type)) { 
    temp<-strsplit(as.character(read_counts$sample[i]), split = 'final') [[1]][1]
    if( temp == combineddf$type[j]) { 
      combineddf$rpkm[j]= (10e9*combineddf$occurences[j]) / (read_counts$count[i] * 3088286401)      
    }
  }
}

cmv_rpkm_plot<-ggplot(combineddf, aes( x=combineddf$length)) + 
  geom_histogram(binwidth = 5, aes()) + 
  #  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")

cmv_rpkm_plot

ggplot(combineddf, aes(y=combineddf$rpkm, x=combineddf$length)) + 
  geom_point() + 
  geom_vline(xintercept = mean(combineddf$length)) +
  xlim(0,500) + 
  ylab('occurrences/total read count of sample')+
  xlab('length')+
  labs(colour="file") + 
  theme_classic() +
  theme(legend.position="none")
######
binned_occurences<-c()
binned_lengths<-c()
for( i in 1:nrow(combineddf)){ 
  if(i%%5==0) {
    tempcount=0
    tempcount<-sum(combineddf$occurences[i],combineddf$occurences[i-1],combineddf$occurences[i-2],combineddf$occurences[i-3],combineddf$occurences[i-4],combineddf$occurences[i-5])
    binned_occurences<-append(binned_occurences,tempcount)
    binned_lengths<-append(binned_lengths,i)
    }
}
#binned_lengths<-as.character(binned_lengths)
#binned_lengths<-as.character(binned_occurences)
binned_df<-cbind.data.frame(binned_occurences,binned_lengths)




######
##########
for(i in 1:nrow(read_counts)) { 
  if(read_counts$cmv_ct[i]==0) { 
    read_counts$cmv_ct[i]=NA
    }
  }
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x,decimals))
}
read_counts$cmvnorm = 100*(read_counts$cmv_quant*250000)/(read_counts$bglobin_quant*3000000000)
percentvspercent<-ggplot(read_counts, aes( x=read_counts$cmvnorm,y= read_counts$percentage, color= read_counts$sample)) + 
  geom_point(size=6) + 
  scale_y_continuous(labels=function(x){sprintf("%.4f", x)}) +
  scale_x_continuous(labels=function(x){sprintf("%.6f", x)}) +
  ylab('Percent of cfDNA reads mapping to CMV ') +
  xlab('expected CMV percentage based on qPCR')+
  labs(colour="file") +
  theme(legend.position="right")
percentvspercent
savePlot(percentvspercent, 'normalized_cmv_and_genome_by_ct.pdf')





cmvctbysum<-ggplot(read_counts, aes( x=read_counts$cmv_ct,y= read_counts$summed, color= read_counts$sample)) + 
  geom_point(size=6) + 
  ylab('number of verified reads')+
  xlab('cmv ct')+
  labs(colour="file") +
  theme(legend.position="right")
cmvctbysum

savePlot(cmvctbysum, 'cmv_by_ct_vs_sum.pdf')

percentcmvbysum<-ggplot(read_counts, aes( x=read_counts$cmv_ct,y= read_counts$percentage, color= read_counts$sample)) + 
  geom_point(size=6) + 
  xlab('cmv ct')+
  ylab('percent of cmv total reads going to cmv ')+
  scale_y_continuous(labels=function(x){sprintf("%.4f", x)}) +
  labs(colour="file") +
  theme(legend.position="right")
percentcmvbysum
savePlot(percentcmvbysum, 'cmv_by_ct_vs_percent.pdf')




#cmv and bglobin quants in copies/ml
qpcr_results<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/qpcr_results.xls', sheetIndex = 3, header = TRUE)

read_counts$cmv_bglobin_ratio_ml<-NA
for(i in 1:nrow(qpcr_results)){ 
  for(j in 1:nrow(read_counts)){
    if(grepl(as.character(qpcr_results$Sample.Name[i]), as.character(read_counts$sample[j]))){ 
      read_counts$cmv_bglobin_ratio_ml[j]<-qpcr_results$CMv_quantity.ml_plasma[i]/qpcr_results$bglobin_quantity.ml_plasma[i]
      
    }
  }
}

figure1b<-ggplot(read_counts, aes(x=read_counts$cmv_bglobin_ratio_ml,y=read_counts$percentage,color = read_counts$sample)) + 
  theme_classic() + 
  #xlim(c(0,.4)) + 
  theme(legend.position="none") +
  xlab('CMV copies/mL Beta-globin copies/mL') + 
  ylab('Percentage of reads mapping to CMV per sample') +
  scale_x_continuous(limits = c(0, .4))+  
  scale_y_log10() + 
  geom_point(size = 3)
figure1b                 
ggsave("figure1b_final.pdf", plot = figure1b, width = 4, height = 4, useDingbats= FALSE)


####resequenced sample####
resequenced<-read.table('/Volumes/Seagate8Tb1/resquenced_121R04_D01_CFFv1_NB0222/read_length_counts.txt')

plot<-ggplot(resequenced, aes( x=resequenced$V2, y= resequenced$V1)) + 
  geom_point() + 
  #  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = median(resequenced$V2)) +
  xlim(0,500) + 
  #ylim(0,200) + 
  ylab('occurences')+
  xlab('length')+
  labs(colour="file") +
  theme_classic() + 
  theme(legend.position="bottom")
plot



resequenced_bam<-scanBam('/Users/gerbix/Documents/vikas/NIPT/21419_download/resequenced_121R04_D01_CFFv1_NB0222.bam')






