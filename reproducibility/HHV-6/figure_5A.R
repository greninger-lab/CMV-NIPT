#make sure to include hhv6 aligned to betaglobin, edar, and rpp30
#NEED TO FILL IN START/STOP FOR EACH GENE
library(Biostrings)
library(Rsamtools)
library(GenomicAlignments)
library(edgeR)
library(rtracklayer)
library(openxlsx)
library(ggplot2)
library(gtools)
library(plotrix)
library(foreach)
library(doParallel)
library(seqinr)
library(scales)
library(wesanderson)
library(RColorBrewer)

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6')

data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/read_counts_all_deduplicated.csv')

edarpath<-'/Users/gerbix/Documents/vikas/NIPT/hhv6/figure_5a/edar/bams/sorted/deduplicated'
betaglobinpath<-'/Users/gerbix/Documents/vikas/NIPT/hhv6/figure_5a/bglobin/bams/sorted/deduplicated'
hhv6apath<-'/Users/gerbix/Documents/vikas/NIPT/hhv6/figure_5a/6a_verified_sams/bams'
hhv6bpath<-'/Users/gerbix/Documents/vikas/NIPT/hhv6/figure_5a/6b_verified_sams/bams'


depthcounter<-function(paths,start,stop) {
  list<-c()
  filenames<-list.files(path=paths, pattern='*.bam$')
  for( i in 1:length(filenames)) {
    count<-0
    tempbam<-scanBam(paste0(paths,'/',filenames[i]))
    for( j in 1:length(tempbam[[1]]$pos)){ 
      if(tempbam[[1]]$pos[j]>start & tempbam[[1]]$pos[j]<stop) { 
        count<-count+1
      }
    }
    print(count)
    list[i]<-count
    print(100*(i/length(filenames)))
  }
  combined<-c()
  combined[[1]]<-list
  combined[[2]]<-filenames
  return(combined)
}
bglobinstart<-1545
bglobinstop<-3871
bglobin<-depthcounter(betaglobinpath,bglobinstart,bglobinstop)
bglobincounts<-bglobin[1]
bglobinnames<-bglobin[2]
bglobindf<-data.frame(bglobincounts,bglobinnames)
bglobindf$type<-'Beta-globin'
colnames(bglobindf)[1]<-'counts'
colnames(bglobindf)[2]<-'names'
bglobindf$depth<-(bglobindf$counts * 76)/(bglobinstop-bglobinstart)
#bglobindf$total_reads<-0
bglobindf$read_counts<-NA
for ( i in 1:nrow(bglobindf)){ 
  tempname<-strsplit(as.character(bglobindf$names[i]), '_bglobin')[[1]][1]
  print(tempname)
  bglobindf$read_counts[i]<-data$count[which(grepl(tempname,data$sample))]
}
bglobindf$length<-bglobinstop-bglobinstart

edarstart<-8390
edarstop<-10900
edar<-depthcounter(edarpath,edarstart,edarstop)
edarcounts<-edar[1]
edarnames<-edar[2]
edardf<-data.frame(edarcounts,edarnames)
edardf$type<-'EDAR'
colnames(edardf)[1]<-'counts'
colnames(edardf)[2]<-'names'
edardf$depth<-(edardf$counts * 76)/(edarstop-edarstart)
edardf$read_counts<-NA
for ( i in 1:nrow(edardf)){ 
  tempname<-strsplit(as.character(edardf$names[i]), '_edar')[[1]][1]
  print(tempname)
  edardf$read_counts[i]<-data$count[which(grepl(tempname,data$sample))]
}
edardf$length<- edarstop - edarstart

hhv6astart<-0
hhv6astop<-159000
hhv6a<-depthcounter(hhv6apath,hhv6astart,hhv6astop)
hhv6acounts<-hhv6a[1]
hhv6anames<-hhv6a[2]
hhv6adf<-data.frame(hhv6acounts,hhv6anames)
hhv6adf$type<-'HHV-6A'
colnames(hhv6adf)[1]<-'counts'
colnames(hhv6adf)[2]<-'names'
hhv6adf$depth<-(hhv6adf$counts * 76)/(hhv6astop-hhv6astart)
hhv6adf$read_counts<-NA
for ( i in 1:nrow(hhv6adf)){ 
  tempname<-strsplit(as.character(hhv6adf$names[i]), '.filtered')[[1]][1]
  print(tempname)
  hhv6adf$read_counts[i]<-data$count[which(grepl(tempname,data$sample))]
}
hhv6adf$length<- hhv6astop - hhv6astart


hhv6bstart<-0
hhv6bstop<-162000
hhv6b<-depthcounter(hhv6bpath,hhv6bstart,hhv6bstop)
hhv6bcounts<-hhv6b[1]
hhv6bnames<-hhv6b[2]
hhv6bdf<-data.frame(hhv6bcounts,hhv6bnames)
hhv6bdf$type<-'HHV-6B'
colnames(hhv6bdf)[1]<-'counts'
colnames(hhv6bdf)[2]<-'names'
hhv6bdf$depth<-(hhv6bdf$counts * 76)/(hhv6bstop-hhv6bstart)
hhv6bdf$read_counts<-NA
for ( i in 1:nrow(hhv6bdf)){ 
  tempname<-strsplit(as.character(hhv6bdf$names[i]), '.filtered')[[1]][1]
  print(tempname)
  hhv6bdf$read_counts[i]<-data$count[which(grepl(tempname,data$sample))]
}
hhv6bdf$length<- hhv6bstop - hhv6bstart


allcombined<-rbind(bglobindf,edardf,hhv6adf,hhv6bdf)

allcombined$normalized<-allcombined$depth/allcombined$read_counts
allcombined$normalized_adjusted<-100*(allcombined$normalized/max(allcombined$normalized))

allcombined$rpkm<-(allcombined$counts * 1e6 * 1e3) / (allcombined$read_counts * (allcombined$length)) 

#not normalized depth
p2 <- ggplot(allcombined, aes(x=factor(type),y=depth, color=allcombined$read_counts))+
  geom_point() + labs(title="Average Coverage") + 
  theme_bw() +
  expand_limits(x = 0, y = 0)

p2

#normalized depth
p3 <- ggplot(allcombined, aes(x=factor(type),y=normalized_adjusted, color=allcombined$type))+
  geom_jitter(width = .3, size = .1) +
  labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  xlab('gene') +
  ylab('normalized depth') +
  scale_y_continuous(trans='log10') +
  theme(legend.position = 'none')
p3

p4 <- ggplot(allcombined, aes(x=factor(type),y=normalized_adjusted, color=allcombined$read_counts))+
  geom_point() +
  labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  xlab('gene') +
  ylab('normalized depth') +
  scale_y_continuous(trans='log10') 
p4

allcombined$fraction<-allcombined$counts/allcombined$read_counts
write.csv(allcombined, 'hhv6_cfdna_info.csv')


#with rpkm 
allcombined$rpkm_adjusted<-allcombined$rpkm / max(allcombined$rpkm)
allcombined$color<-NA
allcombined$color[allcombined$rpkm_adjusted < .1] = 'blue'
allcombined$color[allcombined$rpkm_adjusted > .1  & allcombined$rpkm_adjusted < .8] = 'red'
allcombined$color[allcombined$type=="Beta-globin" | allcombined$type=="EDAR"]<-'green'


p5 <- ggplot(allcombined, aes(x=factor(type),y=rpkm_adjusted, color=allcombined$color))+
  scale_color_manual(values = c( '#436EEE','#077524', '#FF4500')) + 
  geom_jitter(width = .3, size = 1) +
  #labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  theme(legend.position = 'none') + 
  #xlab('') +
  ylab('normalized depth') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),axis.title.x=element_blank())+ 
  scale_y_continuous(trans='log10') +
  theme(text = element_text(size=8)) 
p5
ggsave('figure_5A_full_hhv6_genome.pdf', p5, width = 3, height = 3)



## EDITS

ab_combined<-allcombined
ab_combined$shape<-'circle'
for(i in 1:nrow(ab_combined)){ 
  if(ab_combined$type[i] == 'HHV-6B'){ 
    ab_combined$shape[i]<-'HHV-6B'
  }
  if(ab_combined$type[i] == 'HHV-6A'){ 
    ab_combined$shape[i]<-'HHV-6A'
    }
}
ab_combined$type[ab_combined$type == 'HHV-6A' | ab_combined$type == 'HHV-6B'] <-'HHV-6'

to_remove<-c()
ab_combined$type<-as.character(ab_combined$type)
ab_combined$names<-as.character(ab_combined$names)
hhv6_unique<-unique(ab_combined$names[ab_combined$type == 'HHV-6'])

for(i in 1:length(hhv6_unique)){ 
  temp<-which(ab_combined$names == hhv6_unique[i] & ab_combined$type == 'HHV-6')
  if(ab_combined$rpkm[temp[1]] > ab_combined$rpkm[temp[2]]){ 
    to_remove<-append(to_remove, temp [2])
  }
  else{ 
    to_remove<-append(to_remove, temp[1])
    }
  }
ab_combined<-ab_combined[-to_remove,]


p6 <- ggplot(ab_combined, aes(x=factor(type),y=rpkm_adjusted, color=ab_combined$color))+
  scale_color_manual(values = c( '#436EEE','#077524', '#FF4500')) + 
  geom_jitter(width = .3, size = 1, aes(shape = ab_combined$shape)) +
  #labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  theme(legend.position = 'none') + 
  #xlab('') +
  ylab('normalized depth') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),axis.title.x=element_blank())+ 
  scale_y_continuous(trans='log10') +
  theme(text = element_text(size=8)) 
p6
ggsave('figure_5A_full_hhv6_shape_combined.pdf', p6, width = 3, height = 3)


save.image(file = 'figure_5A.rdata')


