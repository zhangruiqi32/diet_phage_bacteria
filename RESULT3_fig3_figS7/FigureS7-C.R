##Script to Supplementary fig7-C.

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

rongyuan_reads_chayi4 <- aggregate(list(Contig_Abundance=rongyuan_reads_chayi$mapped_Normalized_Mapped_reads) ,
                                   by=list(
                                     PhageLifestyle=rongyuan_reads_chayi$Phacts_result,
                                     mapped_Group=rongyuan_reads_chayi$mapped_Group,
                                     mapped_Contig=rongyuan_reads_chayi$mapped_Contig,
                                     mapped_Diet=rongyuan_reads_chayi$mapped_Diet,
                                     mapped_GD=rongyuan_reads_chayi$mapped_GD),
                                   median )
rongyuan_reads_chayi4 <- spread(rongyuan_reads_chayi4,key =  'PhageLifestyle', value = 'Contig_Abundance')
rongyuan_reads_chayi4[is.na(rongyuan_reads_chayi4)] <- 0
rongyuan_reads_chayi4$Lytic_chu_Temperate <- (as.numeric(rongyuan_reads_chayi4$Lytic)+0.1)/(as.numeric(rongyuan_reads_chayi4$Temperate)+0.2+as.numeric(rongyuan_reads_chayi4$Lytic))
rongyuan_reads_chayi4$Lytic_chu_Temperate <- (as.numeric(rongyuan_reads_chayi4$Lytic)+1)/(as.numeric(rongyuan_reads_chayi4$Temperate)+1)
rongyuan_reads_chayi4$Lytic_chu_Temperate <- log10(rongyuan_reads_chayi4$Lytic_chu_Temperate)
rongyuan_reads_chayi4 <- rongyuan_reads_chayi4[rongyuan_reads_chayi4$mapped_GD%in%c("D0","D10"),]
rongyuan_reads_chayi4$mapped_Group <- factor(rongyuan_reads_chayi4$mapped_Group,
                                             levels = c("CON_D0", "CON_D10","HFD_D0", "HFD_D10","FUC_D0", "FUC_D10"),
                                             labels = c("CON_D0", "CON_D10","HFD_D0", "HFD_D10","FUC_D0", "FUC_D10")  )


############################## [FigS7-C] Painting #####################################
ggplot(rongyuan_reads_chayi4,aes(x=mapped_Group,y=Lytic_chu_Temperate,fill=mapped_Group))+geom_boxplot()+
  geom_signif(  
    comparisons = list( c("FUC_D0", "FUC_D10") ,c("HFD_D0", "HFD_D10") , c("CON_D0", "CON_D10") ) )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("")+scale_fill_manual(values = 
                               c("#DDDDDD","#8BBCD6",
                                 "#DDDDDD","#DD7F88",
                                 "#DDDDDD","#ACD26A"))+
  ylab("Phage lysis and lysogenic ratio")
ggsave("../figS7-C.pdf")



