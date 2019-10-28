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

setwd('/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6')
data <- read.csv("/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/read_counts_all.csv")

#update for deduplication
data<-read.csv('/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/read_counts_all_deduplicated.csv')

# data <- read.xlsx("gtex_icihhv6_hhv6b_positives.xlsx", sheetIndex = 2)
# datatrimmed <- data[, -c(5:9),(12:13)]
# data<-data[complete.cases(data$Run), ]
# data<-data[data$LibrarySelection=='cDNA',]
# #reads z29 gff3, trims out repeats
# gff3<-readGFF('hhv6areference.gff3',version=3,columns=NULL, tags=NULL, filter=NULL, nrows=-1, raw_data=FALSE)

edarpath<-'/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/edar'
rpp30path<-'rpp30'
betaglobinpath<-'/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/bglobin'
hhv6apath<-'/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/hhv6a'
hhv6bpath<-'/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/hhv6b'


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
bglobindf$type<-'Beta globin'
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

hhv6astart<-42000
hhv6astop<-90000
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
  tempname<-strsplit(as.character(hhv6adf$names[i]), '.sam')[[1]][1]
  print(tempname)
  hhv6adf$read_counts[i]<-data$count[which(grepl(tempname,data$sample))]
}
hhv6adf$length<- hhv6astop - hhv6astart


hhv6bstart<-42000
hhv6bstop<-90000
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
  tempname<-strsplit(as.character(hhv6bdf$names[i]), '.sam')[[1]][1]
  print(tempname)
  hhv6bdf$read_counts[i]<-data$count[which(grepl(tempname,data$sample))]
}
hhv6bdf$length<- hhv6bstop - hhv6bstart

#allcombined<-rbind(bglobindf,edardf,rpp30df,hhv6adf,hhv6bdf)
allcombined<-rbind(bglobindf,edardf,hhv6adf,hhv6bdf)
#allcombined$total_reads<-allcombined$total_reads*1000000
allcombined$normalized<-allcombined$depth/allcombined$read_counts
allcombined$normalized_adjusted<-100*(allcombined$normalized/max(allcombined$normalized))
#allcombined$normalized<-allcombined$depth/allcombined$total_reads
#allcombined$normalized<-allcombined$normalized*300000000000
#allcombined[nrow(allcombined)+1,] <- NA
#allcombined[nrow(allcombined),3]<-'negative'
#allcombined[nrow(allcombined),4:6]<-0

allcombined$rpkm<-(allcombined$counts * 1e6 * 1e3) / (allcombined$read_counts * (allcombined$length)) 

pal <- wes.palette(name = "Zissou", type = "continuous")

graph<-ggplot(allcombined,aes(x=names, y=depth, color=as.factor(type))) + geom_point()
graph + scale_color_brewer(palette="Dark2") + theme_minimal() +
  labs(x = "Gene")

#not normalized depth
p2 <- ggplot(allcombined, aes(x=factor(type),y=depth, color=allcombined$read_counts))+
  geom_point() + labs(title="Average Coverage") + 
  theme_bw() +
  expand_limits(x = 0, y = 0)

p2

ggsave('nipt_cfdna_average_coverage.pdf', p2)




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
ggsave('nipt_cfdna_average_coverage_normalized.pdf', height = 3, width = 3)

p4 <- ggplot(allcombined, aes(x=factor(type),y=normalized_adjusted, color=allcombined$read_counts))+
  geom_point() +
  labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  xlab('gene') +
  ylab('normalized depth') +
  scale_y_continuous(trans='log10') 
p4
ggsave('nipt_cfdna_average_coverage_colored.pdf', p4)

allcombined$fraction<-allcombined$counts/allcombined$read_counts
write.csv(allcombined, 'hhv6_cfdna_info.csv')


#with rpkm 
allcombined$rpkm_adjusted<-allcombined$rpkm / max(allcombined$rpkm)
p5 <- ggplot(allcombined, aes(x=factor(type),y=rpkm_adjusted, color=allcombined$type))+
  geom_jitter(width = .3, size = 1) +
  labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  theme(legend.position = 'none') + 
  xlab('gene') +
  ylab('normalized depth') +
  scale_y_continuous(trans='log10') 
p5
ggsave('hhv6_nipt_cfdna_average_coverage_rpkm_deduplicated.pdf', p5, width = 3, height = 3)

  allcombined$sample<-NA
for(i in 1:nrow(allcombined)){
  if(grepl('bglobin',allcombined$names[i])){
    allcombined$sample[i]<-strsplit(as.character(allcombined$names[i]), '_bglobin*')[[1]][1]
  }
  else if (grepl('edar',allcombined$names[i])){
    allcombined$sample[i]<-strsplit(as.character(allcombined$names[i]), '_edar*')[[1]][1]
  }
  else { 
    allcombined$sample[i]<-strsplit(as.character(allcombined$names[i]), '.sam*')[[1]][1]
    
    }
  }
to_change_shape<-c()
for(i in 1:nrow(allcombined)){ 
  if(allcombined$rpkm[i] < .2) { 
    to_change_shape<-append(to_change_shape,allcombined$sample[i])
    }
}
to_change_shape<-unique(to_change_shape)
allcombined$shape<-'hhv6b rpkm > .2'

for(i in 1:nrow(allcombined)){
  for(j in 1:length(to_change_shape)){ 
    if(to_change_shape[j] == allcombined$sample[i]){ 
      allcombined$shape[i]<-'hhv6b rpkm < .2'
      }
    }
  }

p6 <- ggplot(allcombined, aes(x=factor(type),y=rpkm, color=allcombined$type, shape = allcombined$shape))+
  geom_jitter(width = .3) +
  labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  xlab('gene') +
  ylab('normalized depth') +
  scale_y_continuous(trans='log10') 
p6
ggsave('nipt_cfdna_average_coverage_rpkm_by_shape.pdf', p6)


no_hhv6a<-allcombined[allcombined$type!='HHV-6A',]

p7 <- ggplot(no_hhv6a, aes(x=factor(type),y=rpkm))+
  geom_jitter(width = .3, size = .5) +
  labs(title="Average Coverage", color = 'Target') + 
  theme_bw() +
  xlab('gene') +
  ylab('normalized depth') +
  scale_y_continuous(trans='log10') 
p7
ggsave('nipt_cfdna_average_coverage_rpkm.pdf', p5)







#fetal fraction data 
write.csv(just_hhv6b,'just_icihhv-6b.csv')
fetal_fraction<-read.xlsx('/Users/gerbix/Documents/vikas/NIPT/31119_download/hhv6/icihhv6_cellfreefetal_added.xlsx')
just_hhv6b<-allcombined[allcombined$type=='HHV-6B',]
just_hhv6b$sample<-as.character(just_hhv6b$sample)
just_hhv6b$fraction<-NA
just_hhv6b$sample<-as.character(just_hhv6b$sample)
fetal_fraction$X1<-as.character(fetal_fraction$X1)
for(i in 1:nrow(just_hhv6b)){ 
  for(j in 1:nrow(fetal_fraction)){ 
    if(grepl(just_hhv6b$sample[i],fetal_fraction$X1[j])){ 
      just_hhv6b$fraction[i]<-fetal_fraction$Predicted_Fetal_Fraction[j]
      }
    }
  }


fetal_fraction_plot<-ggplot(just_hhv6b, aes(x = just_hhv6b$fraction, y= just_hhv6b$rpkm, group = just_hhv6b$shape)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE) + 
  theme_classic() + 
  ylab('RPKM') + 
  xlab('Predicted fetal fraction')
fetal_fraction_plot
ggsave(plot = fetal_fraction_plot,'icihHV-6_cfdna_fetalfractions.pdf', height = 3, width = 3)

male_fetal_fraction<-fetal_fraction[fetal_fraction$Sex_Classification=='XY',]
just_hhv6b_male<-allcombined[allcombined$type=='HHV-6B',]
just_hhv6b_male$sample<-as.character(just_hhv6b_male$sample)
just_hhv6b_male$fraction<-NA
just_hhv6b_male$sample<-as.character(just_hhv6b_male$sample)
male_fetal_fraction$X1<-as.character(male_fetal_fraction$X1)
for(i in 1:nrow(just_hhv6b_male)){ 
  for(j in 1:nrow(male_fetal_fraction)){ 
    if(grepl(just_hhv6b_male$sample[i],male_fetal_fraction$X1[j])){ 
      just_hhv6b_male$fraction[i]<-male_fetal_fraction$Predicted_Fetal_Fraction[j]
    }
  }
}
male_fetal_fraction_plot<-ggplot(just_hhv6b_male, aes(x = just_hhv6b_male$fraction, y= just_hhv6b_male$rpkm, group = just_hhv6b_male$shape)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE) + 
  theme_classic() + 
  ylab('RPKM') + 
  xlab('Predicted fetal fraction')
male_fetal_fraction_plot
ggsave(plot = male_fetal_fraction_plot,'male_icihHV-6_cfdna_fetalfractions.pdf', height = 3, width = 3)









