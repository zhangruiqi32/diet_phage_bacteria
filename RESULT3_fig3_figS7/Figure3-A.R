##Script to Figure 3-A.

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
virsorter_contents <- read.csv("~/cooperation/202409zhaofengxiang/data/fig3A_virsorter_contents.csv")
checkv_contents <- read.csv("~/cooperation/202409zhaofengxiang/data/fig3A_checkv_contents.csv")
viruse_contig <- read.csv("~/cooperation/202409zhaofengxiang/data/fig3A_viruse.csv")
#viruse_contig <- read.table("/home/wanglab/database/diet/Fucoidan_5.4/viruse_contig2")


############################## Data processing ##############################
viruse_contig$seqname <- sapply(strsplit(viruse_contig$Contig_number,split = "|",fixed = T),"[",1 )
viruse_contig2 <- merge(viruse_contig ,virsorter_contents[c("seqname","group","vir","shape")],by="seqname")
viruse_contig2 <- merge(viruse_contig2 ,checkv_contents,by.x="Contig_number",by.y="contig_id")

viruse_contig3 <-viruse_contig2[is.na(viruse_contig2$viruse_family)==F,]
viruse_contig3$viruse_family[viruse_contig3$viruse_family%in%names(table(viruse_contig3$viruse_family)[table(viruse_contig3$viruse_family)<40])] <- "Others"
color<- c(brewer.pal(8,'Set2'),brewer.pal(7,'Paired'))
color_name <- unique(viruse_contig3$viruse_family)
color <- c("#D3D3D3",brewer.pal(8,'Set2')[3:7])
color <- c("#C0C0C0",brewer.pal(8,'Set1')[1:5])
names(color)<-color_name
size_name=unique(viruse_contig2$shape)
size <-c(2,4)
names(size)   <- size_name
alpha <- c(0.6,0.6,0.6,0.6,0.6,0.6)
alpha <- c(0.4,0.5,0.5,0.5,0.5,0.5)
names(alpha) <- unique(viruse_contig3$viruse_family)

############################## [Fig3-A] Painting #####################################
p1_contig<- ggplot(viruse_contig3,aes(x=log10(length),y=log10(covrage),size=shape,color=viruse_family,alpha=viruse_family))+geom_point(shape=21)+theme_bw()+
  # scale_color_manual(values = color)+
  scale_size_manual(values = size)+
  scale_color_manual(values = color)+
  scale_alpha_manual(values =alpha)+theme_bw()+theme(legend.position = "bottom")
p2_contig<- ggplot(viruse_contig3, aes(x = log10(length))) +
  geom_density(fill="gray",alpha=0.3,adjust=1)+theme_bw()+theme(axis.ticks=element_blank(),axis.text=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank())+
  xlab("")+ylab("")
p3_contig<- ggplot(viruse_contig3, aes(y = log10(covrage))) +
  geom_density(fill="gray",alpha=0.3,adjust=1)+theme_bw()+theme(axis.ticks=element_blank(),axis.text=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank())+
  xlab("")+ylab("")


ggarrange(p2_contig,NULL,p1_contig,p3_contig, heights = c(0.5,2),widths = c(2, 0.5),nrow=2,ncol=2,common.legend = TRUE,legend = "right")
