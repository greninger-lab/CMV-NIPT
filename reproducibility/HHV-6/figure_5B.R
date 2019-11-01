library(svMisc)

setwd('/Users/gerbix/Documents/vikas/NIPT/hhv6/figure_5b/6a_verified_sams/bams')

read_counts <- read.csv("/Users/gerbix/Documents/vikas/NIPT/hhv6/figure_5b/read_counts_all_deduplicated.csv")
system("for i in *.bam ; do echo processing $i ; samtools view -@ 8 $i | awk '{print $9}' | sort | uniq -c > $i.results.txt ; done")
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
      #savePlot(plot, plotname)
      dflist[[i]]<-histo
    }
  }
}

setwd('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6')


combined <- do.call("rbind", dflist)
lengthlist<-c()
samplelist<-c()
for(i in 1:nrow(combined)){
  for(j in 1:combined$occurences[i]){
    lengthlist<-append(lengthlist,combined$length[i])
    samplelist<-append(samplelist,combined$type[i])
  }
}


combined_expanded<-data.frame(lengthlist,samplelist)

for(i in 1:nrow(combined_expanded)){ 
  progress(i,nrow(combined_expanded))
combined_expanded$sample_trimmed<-strsplit(as.character(combined_expanded$samplelist),'_')[[i]][1]
}

low_cluster<-c('99P03_C01','120R22_F03','104P04_D01','98P11','97P10_B02','93R20_D03')

combined_expanded$cluster<-"High HHV-6"
for(i in 1:nrow(combined_expanded)){ 
  if(grepl(paste(low_cluster,collapse="|"),combined_expanded$samplelist[i])){ 
    combined_expanded$cluster[i]<-"Low HHV-6"
    }
  }

plot<-ggplot(combined_expanded, aes( x=combined_expanded$lengthlist)) + 
  geom_histogram(binwidth = 5, aes(fill = combined_expanded$cluster, y=(100 * ..count../sum(..count..)))) + 
  scale_fill_manual(values=c("orangered", "royalblue2")) + 
  xlim(0,500) + 
  ylab('Percent within each alignment')+
  xlab('Fragment size')+
  labs(colour="file") +
  theme_classic() +
  #theme(legend.position = 'bottom') +
  labs(fill="") +
  theme(text = element_text(size=10)) +
  theme(legend.position = c(0.8, 0.3)) + 
  theme(legend.title=element_blank())
plot
ggsave(plot=plot, 'figure_5B.pdf', height = 3, width = 3)

write.csv(combined_expanded, 'hhv6_isizes_all.csv')

save.image(file = 'figure_5B.rdata')


