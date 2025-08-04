##Script to Figure 4-BC and Supplementary fig8-A.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","ggsignif","reshape2","ggpubr","ggTimeSeries","RColorBrewer","ggplot2")
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
HGT_reads <- read.csv("../data/fig4BC_S8A_reads.csv")
HGT_contig_length <- read.csv("../data/fig4BC_S8A_length.csv")

viruse_contig <- read.table("../data/viruse_contig")
x <- read.csv("../data/DABs_lefse.csv")
meta <- read.table("../data/sample_meta")

############################## Data processing ##############################
HGT_contig_length$Contig2 <- sapply(strsplit(HGT_contig_length$Contig,split ="-"),"[",1)
HGT_reads<- merge(HGT_reads,HGT_contig_length,by.x="Contig",by.y="Contig2")
HGT_reads$Normalized_Mapped_reads <- as.numeric(HGT_reads$reads_number)/(as.numeric(HGT_reads$HGT_contig_length)*as.numeric(HGT_reads$num_seqs))*10000000000
HGT_reads2 <-spread(HGT_reads[colnames(HGT_reads)%in%c("num_seqs","reads_number","GD","Diet","Group")==F],key=Sample_id,value =Normalized_Mapped_reads )
HGT_reads2[is.na(HGT_reads2)] <- 0
HGT_reads2 <-gather(HGT_reads2,key = "Sample_id",value = "Normalized_Mapped_reads",L1EGH020932:L1HGH010566)
HGT_reads2 <-merge(HGT_reads2,meta,by="Sample_id")


HGT_contigs <- paste(files,"VIBRANT_out2/VIBRANT_viruses_cd-hit2/VIBRANT_phages_viruses_cd-hit2/viruses_cd-hit2.phages_combined.faa",sep="/")
HGT_contig_contents <- data.frame(matrix(ncol=1,nrow=0))
for (HGT_index in 1:length(HGT_contigs)){
  HGT_contig_content<- read.csv(HGT_contigs[HGT_index],sep="\t",header = F)
  if(nrow(HGT_contig_content)!=0){
    HGT_contig_content$Sample_id <- files[HGT_index]
  }
  HGT_contig_contents <- rbind(HGT_contig_contents,HGT_contig_content )
}
HGT_contig_contents <- HGT_contig_contents[grep(">",HGT_contig_contents$V1,fixed = T),]
HGT_contig_contents <- HGT_contig_contents[grep("integrase|recombinase|transposase",HGT_contig_contents$V5,ignore.case = T),]

HGT_contig_contents2 <- HGT_contig_contents[grep("integrase|recombinase|transposase",HGT_contig_contents$V5,ignore.case = T),]
HGT_contig_contents2$scaffold2 <- substring(HGT_contig_contents2$V1, 2, last = 31)
HGT_contig_contents2$type <- apply(HGT_contig_contents2,1,function(x) if ( length(grep("recombinase",x[5],ignore.case = T))!=0) {"recombinase"
} else if (  length(grep("integrase",x[5],ignore.case = T))!=0 ) {
  "integrase"
} else if (  length(grep("transposase",x[5],ignore.case = T))!=0 ) {
  "transposase"
})
HGT_contig_contents2$type <- apply(HGT_contig_contents2,1,function(x) if ( length(grep("integrase",x[5],ignore.case = T))!=0) {"integrase"
} else if (  length(grep("recombinase",x[5],ignore.case = T))!=0 ) {
  "recombinase"
} else if (  length(grep("transposase",x[5],ignore.case = T))!=0 ) {
  "transposase"
})


chayijun <- sapply(strsplit(x$Spe,split = "g__"),"[",2)
chayijun <- paste(chayijun,collapse = "|")
viruse_IMG_contig_anno_chayijun <- viruse_contig[grep(chayijun,viruse_contig$Bacteria_Host),]

viruse_IMG_contig_anno_chayijun$Contig_name2 <- viruse_IMG_contig_anno_chayijun$Contig_number
viruse_IMG_contig_anno_chayijun$Contig_name2 <- substring(viruse_IMG_contig_anno_chayijun$Contig_name2, 1, last = 30)
HGT_viruse_IMG_contig_anno_chayijun <- data.frame(viruse_IMG_contig_anno_chayijun[viruse_IMG_contig_anno_chayijun$Contig_name2%in%HGT_contig_contents2$scaffold2,])


HGT_reads$Contig2 <-  substring(HGT_reads$Contig, 1, last = 30)
HGT_reads3 <- HGT_reads
colnames(HGT_reads3)[which(colnames(HGT_reads3)=="Sample_id")] <- "Mapped_Sample_id"
colnames(HGT_reads3)[which(colnames(HGT_reads3)=="GD")] <- "Mapped_GD"
colnames(HGT_reads3)[which(colnames(HGT_reads3)=="Diet")] <- "Mapped_Diet"
colnames(HGT_reads3)[which(colnames(HGT_reads3)=="Normalized_Mapped_reads")] <- "HGT_Normalized_Mapped_reads"
HGT_viruse_IMG_contig_anno_chayijun <-merge(HGT_viruse_IMG_contig_anno_chayijun,HGT_reads3[c("Contig2","Mapped_Sample_id","HGT_Normalized_Mapped_reads","Mapped_GD","Mapped_Diet")],by.x="Contig_name2",by.y="Contig2")
HGT_viruse_IMG_contig_anno_chayijun$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun$Bacteria_Host,split = "s__"),"[",1)
HGT_viruse_IMG_contig_anno_chayijun$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun$Bacteria_Host_G,split = "g__"),"[",2)
HGT_viruse_IMG_contig_anno_chayijun$fill <- paste(HGT_viruse_IMG_contig_anno_chayijun$Contig_name2,HGT_viruse_IMG_contig_anno_chayijun$Mapped_GD,sep="_")
HGT_viruse_IMG_contig_anno_chayijun$fill <- paste(HGT_viruse_IMG_contig_anno_chayijun$Bacteria_Host_G,HGT_viruse_IMG_contig_anno_chayijun$Mapped_GD,sep="_")
HGT_viruse_IMG_contig_anno_chayijun$Mapped_Group <- paste(HGT_viruse_IMG_contig_anno_chayijun$Mapped_Diet,HGT_viruse_IMG_contig_anno_chayijun$Mapped_GD,sep="_")

HGT_viruse_IMG_contig_anno_chayijun <- HGT_viruse_IMG_contig_anno_chayijun[order(HGT_viruse_IMG_contig_anno_chayijun$Bacteria_Host_G,HGT_viruse_IMG_contig_anno_chayijun$HGT_Normalized_Mapped_reads),]
HGT_viruse_IMG_contig_anno_chayijun$Contig_number <- factor(HGT_viruse_IMG_contig_anno_chayijun$Contig_number,levels = unique(HGT_viruse_IMG_contig_anno_chayijun$Contig_number),
                                                            labels = unique(HGT_viruse_IMG_contig_anno_chayijun$Contig_number) )

HGT_viruse_IMG_contig_anno_chayijun2<- HGT_viruse_IMG_contig_anno_chayijun
HGT_viruse_IMG_contig_anno_chayijun2 <- aggregate(list(mapped_Normalized_Mapped_reads=HGT_viruse_IMG_contig_anno_chayijun2$HGT_Normalized_Mapped_reads),
                                                  by=list(mapped_Group=HGT_viruse_IMG_contig_anno_chayijun2$Mapped_Group,
                                                          Contig_number=HGT_viruse_IMG_contig_anno_chayijun2$Contig_number,
                                                          Bacteria_Host=HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host),median) 
HGT_viruse_IMG_contig_anno_chayijun2$Diet <- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun2$mapped_Group,split = "_"),"[",1)
HGT_viruse_IMG_contig_anno_chayijun2$Time<- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun2$mapped_Group,split = "_"),"[",2)
HGT_viruse_IMG_contig_anno_chayijun2 <- spread(HGT_viruse_IMG_contig_anno_chayijun2[,-c(which(colnames(HGT_viruse_IMG_contig_anno_chayijun2)=="mapped_Group"))], key= 'Time', value =  'mapped_Normalized_Mapped_reads')
HGT_viruse_IMG_contig_anno_chayijun2[is.na(HGT_viruse_IMG_contig_anno_chayijun2)] <- 0
HGT_viruse_IMG_contig_anno_chayijun2$Change <- as.numeric(HGT_viruse_IMG_contig_anno_chayijun2$D10-HGT_viruse_IMG_contig_anno_chayijun2$D0)
HGT_viruse_IMG_contig_anno_chayijun2$Contig_number2 <- substring(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number,1,30)
HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host,split = "s__"),"[",1)
HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G,split = "g__"),"[",2)
HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G,split = ";"),"[",1)
HGT_viruse_IMG_contig_anno_chayijun2 <-HGT_viruse_IMG_contig_anno_chayijun2[order(HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G,HGT_viruse_IMG_contig_anno_chayijun2$Diet,HGT_viruse_IMG_contig_anno_chayijun2$Change),]
HGT_viruse_IMG_contig_anno_chayijun2$Contig_number <- factor(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number,levels = unique(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number),labels = unique(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number))
HGT_viruse_IMG_contig_anno_chayijun2$Change2[as.numeric(HGT_viruse_IMG_contig_anno_chayijun2$Change)<0] <- -1*(log10(abs(HGT_viruse_IMG_contig_anno_chayijun2$Change[as.numeric(HGT_viruse_IMG_contig_anno_chayijun2$Change)<0])))
HGT_viruse_IMG_contig_anno_chayijun2$Change2[as.numeric(HGT_viruse_IMG_contig_anno_chayijun2$Change)>0] <- log10(abs(HGT_viruse_IMG_contig_anno_chayijun2$Change[as.numeric(HGT_viruse_IMG_contig_anno_chayijun2$Change)>0]))
a<- unique(HGT_contig_contents2[c("scaffold2","type")])
a <- a[!duplicated(a$scaffold2),]
HGT_viruse_IMG_contig_anno_chayijun2 <- merge(HGT_viruse_IMG_contig_anno_chayijun2,a,by.x="Contig_number2",by.y="scaffold2")
HGT_viruse_IMG_contig_anno_chayijun2 <- spread(HGT_viruse_IMG_contig_anno_chayijun2[,colnames(HGT_viruse_IMG_contig_anno_chayijun2)%in%c("D0","D10","Change","D14","D18")==F],key = "Diet",value = "Change2")
HGT_viruse_IMG_contig_anno_chayijun2[is.na(HGT_viruse_IMG_contig_anno_chayijun2)] <- 0

HGT_viruse_IMG_contig_anno_chayijun2 <- melt(HGT_viruse_IMG_contig_anno_chayijun2,id.vars = colnames(HGT_viruse_IMG_contig_anno_chayijun2)[1:5],variable.name = "Diet",value.name  = "Change2")
HGT_viruse_IMG_contig_anno_chayijun2$Diet <- factor(HGT_viruse_IMG_contig_anno_chayijun2$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC"))
HGT_viruse_IMG_contig_anno_chayijun2$type2 <- paste(HGT_viruse_IMG_contig_anno_chayijun2$type,factor(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number,levels =unique(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number) ,labels = c(1:length(unique(HGT_viruse_IMG_contig_anno_chayijun2$Contig_number)))),sep = "_")
HGT_viruse_IMG_contig_anno_chayijun2 <-HGT_viruse_IMG_contig_anno_chayijun2[order(HGT_viruse_IMG_contig_anno_chayijun2$type,HGT_viruse_IMG_contig_anno_chayijun2$Bacteria_Host_G,HGT_viruse_IMG_contig_anno_chayijun2$Diet,HGT_viruse_IMG_contig_anno_chayijun2$Change2),]
HGT_viruse_IMG_contig_anno_chayijun2$type2 <- factor(HGT_viruse_IMG_contig_anno_chayijun2$type2,levels = unique(HGT_viruse_IMG_contig_anno_chayijun2$type2),labels = unique(HGT_viruse_IMG_contig_anno_chayijun2$type2))
HGT_viruse_IMG_contig_anno_chayijun2$type2 <- factor(HGT_viruse_IMG_contig_anno_chayijun2$type2,levels = unique(HGT_viruse_IMG_contig_anno_chayijun2$type2),labels = length(unique(HGT_viruse_IMG_contig_anno_chayijun2$type2)):1   )
HGT_viruse_IMG_contig_anno_chayijun2$type2 <- paste0("Contig",HGT_viruse_IMG_contig_anno_chayijun2$type2)
HGT_viruse_IMG_contig_anno_chayijun2$type2<- factor(HGT_viruse_IMG_contig_anno_chayijun2$type2,levels = unique(HGT_viruse_IMG_contig_anno_chayijun2$type2),labels = unique(HGT_viruse_IMG_contig_anno_chayijun2$type2)   )
colnames(HGT_viruse_IMG_contig_anno_chayijun2)[which(colnames(HGT_viruse_IMG_contig_anno_chayijun2)=="Bacteria_Host_G")] <- "Phage"
HGT_viruse_IMG_contig_anno_chayijun2$Phage <- paste(sapply(strsplit(HGT_viruse_IMG_contig_anno_chayijun2$Phage,split=";"),"[",1),"phage",sep = "_")

HGT_viruse_IMG_contig_anno2<-merge(HGT_contig_contents2[c("scaffold2","type")],HGT_reads3[c("Contig2","Mapped_Sample_id","HGT_Normalized_Mapped_reads","Mapped_GD","Mapped_Diet","Contig")],by.x="scaffold2",by.y="Contig2")
HGT_viruse_IMG_contig_anno2<-merge(HGT_viruse_IMG_contig_anno2,viruse_contig,by.x="scaffold2",by.y="Contig_number2")
HGT_viruse_IMG_contig_anno2$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno2$Bacteria_Host,split = "s__"),"[",1)
HGT_viruse_IMG_contig_anno2$Bacteria_Host_G <- sapply(strsplit(HGT_viruse_IMG_contig_anno2$Bacteria_Host_G,split = "g__"),"[",2)

a<- aggregate(list(HGT_Normalized_Mapped_reads=HGT_viruse_IMG_contig_anno2$HGT_Normalized_Mapped_reads),
              by=list(Mapped_GD=HGT_viruse_IMG_contig_anno2$Mapped_GD,
                      Mapped_Diet=HGT_viruse_IMG_contig_anno2$Mapped_Diet,
                      type=HGT_viruse_IMG_contig_anno2$type) ,sum
)


############################## [Fig4-B] Painting #####################################
color<- c("#98BEE5","#BCD5EF","#E1ECF8")
a$Mapped_Diet <- factor(a$Mapped_Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC")  )

ggplot(a,aes(x=Mapped_GD,y=HGT_Normalized_Mapped_reads,group=type,fill=type))+
  geom_area()+
  facet_wrap(.~Mapped_Diet,)+
  scale_fill_manual(values = color)+
  theme_bw()+ theme(panel.grid=element_blank())+
  xlab("Time")+ylab("The Abundance of HGT Related Genes")
ggsave("../fig4-B.pdf")

############################## [Fig4-C] Painting #####################################
p1<- ggplot(HGT_viruse_IMG_contig_anno_chayijun2,aes(x=Change2,y=type2,fill=Phage))+
  geom_bar(stat = "identity")+facet_wrap(.~Diet,nrow = 1)+
  theme_bw()+ theme(panel.grid=element_blank())+
  scale_fill_manual(values = c("#8DC7B8","#E0C092","#8A9CC4","#DA87B6","#F08961"))+
  #theme(axis.ticks.y =element_blank(),axis.text.y =element_blank())+
  theme(axis.ticks.y =element_blank())+
  xlab("Changes in the Abundance of HGT Related Genes")+ylab("")

p2<- ggplot(HGT_viruse_IMG_contig_anno_chayijun2,aes(x=1,y=type2,fill=type))+
  geom_bar(stat="identity",position="dodge",width=1)+
  scale_fill_manual( values = c("#98BEE5","#BCD5EF","#E1ECF8"))+
  theme_bw()+ theme(panel.grid=element_blank())+
  # scale_color_manual(values = annotation_colors3)+
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text=element_blank(),
    legend.position="left"
  )+xlab("")+ylab("")

ggarrange(p2,p1,nrow = 1,widths = c(0.5, 2))
ggsave("../fig4-C.pdf")

############################## [FigS8-A] Processing and Painting #####################################
HGT_viruse_IMG_contig_anno2$chayi <- "Normal bacteria"
HGT_viruse_IMG_contig_anno2$chayi[grep(chayijun,HGT_viruse_IMG_contig_anno2$Bacteria_Host)] <- "Differential bacteria"
a<- HGT_viruse_IMG_contig_anno2[HGT_viruse_IMG_contig_anno2$chayi=="Differential bacteria",]
a$fill2 <- paste(a$Bacteria_Host_G,a$Mapped_Diet,sep="_")
a<- aggregate(list(HGT_Normalized_Mapped_reads=a$HGT_Normalized_Mapped_reads),
              by=list(Mapped_GD=a$Mapped_GD,
                      Mapped_Diet=a$Mapped_Diet,
                      Bacteria_Host_G=a$Bacteria_Host_G) ,sum
)
a$fill <- paste(a$Mapped_Diet,a$Bacteria_Host_G,sep="_")
a1 <- a[a$Mapped_GD=="D0",]
a$HGT_D0 <- factor(a$fill,levels = a1$fill,labels = a1$HGT_Normalized_Mapped_reads)
a$HGT_change <- a$HGT_Normalized_Mapped_reads -  round(as.numeric(as.character(a$HGT_D0)),8)
a$HGT_change[a$HGT_change>0] <- log10(a$HGT_change[a$HGT_change>0]+1)
a$HGT_change[a$HGT_change<0] <- log10(abs(a$HGT_change[a$HGT_change<0])+1)*(-1)

ggplot(a,aes(x=Mapped_GD,y=Bacteria_Host_G,fill=HGT_change))+
  geom_point(shape=21,size=4)+
  theme_test()+
  facet_wrap(.~Mapped_Diet)+
  scale_fill_gradient2(low = "#2C7FB8",mid = "#F7F8F8",high = "#D7301F")
ggsave("../figS8-A.pdf")

