##Script to Supplementary fig7-E.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","RColorBrewer","ggpubr","ggplot2")
bioconductor_packages=c()

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if (!require("BiocManager", quietly = T, warn.conflicts = F)) install.packages("BiocManager")

for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}


############################## Data loading ##############################
viruse_contig <- read.table("~/cooperation/202409zhaofengxiang/data/viruse_contig")
rongyuan_reads <- read.csv("~/cooperation/202409zhaofengxiang/data/fig3_S7_reads.csv")
rongyuan_contig_length <- read.csv("~/cooperation/202409zhaofengxiang/data/fig3_S7_length.csv")
viruse_IMG_contig_anno_chayijun <- read.csv("~/cooperation/202409zhaofengxiang/data/fig3_S7_anno.csv")


############################## Data processing ##############################
rongyuan_reads$Contig2 <- substring(rongyuan_reads$Contig,1,30)
colnames(rongyuan_contig_length) <- c("Contig","rongyuan_contig_length")
rongyuan_contig_length$Contig2 <- sapply(strsplit(rongyuan_contig_length$Contig,split ="-"),"[",1)
rongyuan_reads<- merge(rongyuan_reads,rongyuan_contig_length,by.x="Contig",by.y="Contig2")
rongyuan_reads$Normalized_Mapped_reads <- as.numeric(rongyuan_reads$reads_number)/(as.numeric(rongyuan_reads$rongyuan_contig_length)*as.numeric(rongyuan_reads$num_seqs))*10000000000
colnames(rongyuan_reads) <- paste("mapped",colnames(rongyuan_reads),sep="_")
rongyuan_reads_chayi <- merge(rongyuan_reads,viruse_IMG_contig_anno_chayijun,by.x="mapped_Contig2",by.y="Contig_number2")
rongyuan_reads_chayi$fill <- paste(rongyuan_reads_chayi$Bacteria_Host_G,rongyuan_reads_chayi$GD_mapped,sep="_")


a <- aggregate(list(Contig_Abundance=rongyuan_reads_chayi$mapped_Normalized_Mapped_reads),
               by=list(Time=rongyuan_reads_chayi$mapped_GD,
                       Phacts_result=rongyuan_reads_chayi$Phacts_result,
                       mapped_Contig2=rongyuan_reads_chayi$mapped_Contig2,
                       Diet=rongyuan_reads_chayi$mapped_Diet,
                       Bacteria_Host_G=rongyuan_reads_chayi$Bacteria_Host_G),median)
a <- a[a$Time=="D10",]
a_phact <- aggregate(list(Contig_Abundance=a$Contig_Abundance),
                     by=list(Time=a$Time,
                             Phacts_result=a$Phacts_result,
                             Diet=a$Diet,
                             Bacteria_Host_G=a$Bacteria_Host_G),sum)
a_phact <- spread(a_phact ,key = "Phacts_result",value = "Contig_Abundance")
a_phact <- a_phact[a_phact$Diet=="CON",]
a_phact$bili <- a_phact$Temperate/(a_phact$Temperate+a_phact$Lytic)
a_phact <- a_phact[order(a_phact$bili),]
a$Bacteria_Host_G <- factor(a$Bacteria_Host_G,levels = a_phact$Bacteria_Host_G,labels = a_phact$Bacteria_Host_G)
a$Diet <- factor(a$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC")  )

############################## [FigS7-E] Painting #####################################
color <- c("#CD776D","#DCDDDD")
names(color) <- c("Lytic","Temperate")

ggplot(a,aes(x=Contig_Abundance,y=Bacteria_Host_G,fill=Phacts_result))+
  geom_bar(position = "fill",stat = 'identity')+facet_wrap(.~Diet,)+
  scale_fill_manual(values = color  )+xlab("")+ylab("")+theme_bw()
ggsave("../figS7-E.pdf")