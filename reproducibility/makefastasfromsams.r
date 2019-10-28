# Arg 1 is the path to the folder containing the deduplicated BAMs
args = commandArgs(trailingOnly=TRUE)

library(Rsamtools)
library(seqinr)
library(svMisc)

filenames <- list.files(args[1], pattern ='.bam$')
setwd(args[1])
sequence= c()
identifier= c()
header=c()
for (i in 1:length(filenames)){ 
  lengthfilenames= length(filenames)
  tempbam=scanBam(filenames[i])
  if(length(tempbam[[1]]$seq)>0){
  for(j in 1:length(tempbam[[1]]$seq)){ 
       sequence =  append(sequence,as.character(tempbam[[1]]$seq[j]))
       tempfilename = paste0(filenames[i],'.',j)
       header=append(header, as.character(tempbam[[1]]$qname[j]))
       identifier = append(identifier,tempfilename)
  }
  }
  }
fastadf <- data.frame(identifier, sequence, header)
fastadf$header<-as.character(fastadf$header)
fastadf$identifier<-as.character(fastadf$identifier)


# repeatmasker limits header names to 50 characters so this deletes the non-unique regions of the header
for(i in 1:nrow(fastadf)){
  progress(i, nrow(fastadf))
  if(grepl('NB502000', fastadf$header[i])) {
    fastadf$header[i]<-gsub(pattern = 'NB502000', replacement = '', fastadf$header[i])
    }
  if(grepl('NS500359', fastadf$header[i])) {
    fastadf$header[i]<-gsub(pattern = 'NS500359', replacement = '', fastadf$header[i])
  }
  fastadf$identifier[i]<-gsub('.sam.bam','',fastadf$identifier[i])
  fastadf$header[i]<-substr(fastadf$header[i], 17, stop = nchar(fastadf$header[i]))
}

# writes fastas
x<-c()
print('making fastas')
for ( i in 1:nrow(fastadf)){
  lengthvar=nrow(fastadf)
  tempfilename = paste0(fastadf$identifier[i],'.fasta')
  write.fasta(sequences= fastadf$sequence[i], names = paste0(tempfilename,'-',fastadf$header[i]), file.out= tempfilename)
  if(nchar(paste0(tempfilename,'-',fastadf$header[i]))>50){ 
    x<-append(x,i)
    }
  } 





