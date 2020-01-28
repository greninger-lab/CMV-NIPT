library(ggplot2)
library(grid)
library(nplr)
library(scales)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggpubr")
library("ggpubr")

setwd("~/Documents/Manuscripts/CMV_cfDNA/newversion/simulation/")
Samplefrag <- read.csv("Sample_frag.csv",stringsAsFactors=F)
Multifrag <- read.csv("Multi_frag.csv",stringsAsFactors=F)

playfrag <- gather(Samplefrag,"X50","X100","X150","X200","X250","X300","X350",key="amplicon",value="counts")
playfrag[playfrag=="X50"] <- 50
playfrag[playfrag=="X100"] <- 100
playfrag[playfrag=="X150"] <- 150
playfrag[playfrag=="X200"] <- 200
playfrag[playfrag=="X250"] <- 250
playfrag[playfrag=="X300"] <- 300
playfrag[playfrag=="X350"] <- 350

playfragMulti <- gather(Multifrag,"X50","X100","X150","X200","X250","X300","X350",key="amplicon",value="counts")
playfragMulti[playfragMulti=="X50"] <- 50
playfragMulti[playfragMulti=="X100"] <- 100
playfragMulti[playfragMulti=="X150"] <- 150
playfragMulti[playfragMulti=="X200"] <- 200
playfragMulti[playfragMulti=="X250"] <- 250
playfragMulti[playfragMulti=="X300"] <- 300
playfragMulti[playfragMulti=="X350"] <- 350

playfrag$amplicon <- as.factor(playfrag$amplicon)
playfrag$amplicon <- factor(playfrag$amplicon, levels = c(50,100,150,200,250,300,350))
playfragMulti$amplicon <- factor(playfragMulti$amplicon, levels = c(50,100,150,200,250,300,350))


Sample_plot <- ggplot(playfrag, aes(x=amplicon, y=counts),show.legend = FALSE) + geom_violin(fill="dark blue") + 
  xlab(label="Amplicon length") + ylab(label="Measured qPCR copies") + theme_classic() + 
  scale_y_continuous(breaks=c(0,100,200,300,400,500),limits=c(0, 560))

Multi_plot <- ggplot(playfragMulti, aes(x=amplicon, y=counts),show.legend = FALSE) + geom_violin(fill="dark green") + 
  xlab(label="Amplicon length") + ylab(label="Measured qPCR copies") + theme_classic() + 
  scale_y_continuous(breaks=c(0,100,200,300,400,500),limits=c(0, 560)) 

figure <- ggarrange(Sample_plot, Multi_plot,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure
ggsave("Figure5_violinplot.pdf",plot=figure,width=6.2,height=3)     


Multi_plot <- ggplot(playfragMulti, aes(x=amplicon, y=counts),show.legend = FALSE) + geom_violin(fill="dark green") + 
  xlab(label="Amplicon length") + ylab(label="Measured qPCR copies") + theme_classic() + #scale_y_log10()
  scale_y_log10(breaks=c(10,30,100,300,1000),limits=c(5, 1000)) 

ggplot(playfragMulti, aes(x=amplicon, y=counts),show.legend = FALSE) + geom_violin(fill="dark green") + 
  xlab(label="Amplicon length") + ylab(label="Measured qPCR copies") + theme_classic() + theme(aspect.ratio=1) +
  scale_y_log10(breaks=c(3,10,30,100,300,1000,3000,10000),labels=c("-2.00","-1.50","-1.00","-0.50","0.00","0.50","1.00","1.50"),limits=c(1, 10000))



scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1",
                          "2" = "Dose 2")

ggsave("logplot.pdf",width=3,height=3)     

