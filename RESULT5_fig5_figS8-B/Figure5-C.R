##Script to Figure 5-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("RColorBrewer","ggplot2")
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
AMG_contents_metabolism2 <- read.csv("../data/fig5C_metabolism.csv")
AMG_contig_length <- read.csv("../data/fig5_contig_length.csv")
AMG_reads <- read.csv("../data/fig5_AMG_reads.csv")
viruse_IMG_contig_anno_chayijun <- read.csv("../data/fig5C_contig_anno_DiffBact.csv")

############################## Data processing ##############################
AMG_reads_chayi  <- AMG_reads[grep("D0|D10",AMG_reads$mapped_GD),]
AMG_reads_chayi2 <- aggregate(list(mapped_Normalized_Mapped_reads=AMG_reads_chayi$mapped_Normalized_Mapped_reads),by=list(mapped_Group=AMG_reads_chayi$mapped_Group,Contig_number=AMG_reads_chayi$mapped_Contig),median) 
AMG_reads_chayi2$Diet <- sapply(strsplit(AMG_reads_chayi2$mapped_Group,split = "_"),"[",1)
AMG_reads_chayi2$Time<- sapply(strsplit(AMG_reads_chayi2$mapped_Group,split = "_"),"[",2)
AMG_reads_chayi2$Contig_number2 <- substring(AMG_reads_chayi2$Contig_number,1,30)
viruse_IMG_contig_anno_chayijun$Bacteria_Host_G <- sapply(strsplit(viruse_IMG_contig_anno_chayijun$Bacteria_Host,split = ";g__"),"[",2)
viruse_IMG_contig_anno_chayijun$Bacteria_Host_G <- sapply(strsplit(viruse_IMG_contig_anno_chayijun$Bacteria_Host_G,split = ";s__"),"[",1)
AMG_reads_chayi2 <-merge(AMG_reads_chayi2,viruse_IMG_contig_anno_chayijun[c("Contig_number2","Bacteria_Host_G")],by.x="Contig_number2",by.y="Contig_number2")
# AMG_contents_metabolism2 <- AMG_contents_metabolism
colnames(AMG_contents_metabolism2) <- paste(colnames(AMG_contents_metabolism2),"source",sep="_")
# AMG_contents_metabolism2[c("scaffold_source","metabolism_source")]
AMG_contents_metabolism2<- unique(AMG_contents_metabolism2[c("protein_source","metabolism_source")])
# AMG_contents_metabolism2<- AMG_contents_metabolism2[!duplicated(AMG_contents_metabolism2$protein_source),]
AMG_reads_chayi2 <- merge(AMG_reads_chayi2,AMG_contents_metabolism2,by.x="Contig_number",by.y="protein_source")
AMG_reads_chayi2$AMG_anno <- paste(1:nrow(AMG_reads_chayi2),AMG_reads_chayi2$metabolism_source,sep="_")
AMG_reads_chayi2 <- AMG_reads_chayi2[order(AMG_reads_chayi2$metabolism_source,AMG_reads_chayi2$Contig_number,AMG_reads_chayi2$Diet),]
# aggregate(AMG_reads_chayi2$Contig_number,by=list(metabolism_source=AMG_reads_chayi2$metabolism_source,AMG_reads_chayi2$Diet),function(x) 1:length(x))
Contig_metabolism<- unique(AMG_reads_chayi2[c("Contig_number","metabolism_source")])
Contig_metabolism$metabolism_name <- paste(Contig_metabolism$metabolism_source,unlist(tapply(Contig_metabolism$Contig_number,Contig_metabolism$metabolism_source,function(x){return(c(1:length(x)))} )) ,sep="_")
AMG_reads_chayi2$AMG_anno <- factor(AMG_reads_chayi2$Contig_number,levels = Contig_metabolism$Contig_number,labels = Contig_metabolism$metabolism_name)
AMG_reads_chayi2 <-AMG_reads_chayi2[order(AMG_reads_chayi2$Bacteria_Host_G,AMG_reads_chayi2$AMG_anno,AMG_reads_chayi2$Diet,AMG_reads_chayi2$mapped_Normalized_Mapped_reads),]
AMG_reads_chayi2$AMG_anno <- factor(AMG_reads_chayi2$AMG_anno,levels = unique(AMG_reads_chayi2$AMG_anno),labels = unique(AMG_reads_chayi2$AMG_anno))
AMG_reads_chayi2$Bacteria_Host_G <- paste(sapply(strsplit(AMG_reads_chayi2$Bacteria_Host_G,split = ";"),"[",1),"phage")


AMG_reads_chayi3 <- AMG_reads_chayi2
AMG_reads_chayi3 <-aggregate(list(mapped_Normalized_Mapped_reads=AMG_reads_chayi3$mapped_Normalized_Mapped_reads),by=list(mapped_Group=AMG_reads_chayi3$mapped_Group,metabolism_source=AMG_reads_chayi3$metabolism_source,Bacteria_Host_G=AMG_reads_chayi3$Bacteria_Host_G),sum)
AMG_reads_chayi3$Group2 <- paste(AMG_reads_chayi3$mapped_Group,AMG_reads_chayi3$metabolism_source,sep="_")
AMG_reads_chayi3 <- AMG_reads_chayi3[order(AMG_reads_chayi3$Group2,AMG_reads_chayi3$Bacteria_Host_G),]
i=1
AMG_reads_chayi_sum <- aggregate(list(mapped_Normalized_Mapped_reads=AMG_reads_chayi3$mapped_Normalized_Mapped_reads),by=list(AMG_reads_chayi3$metabolism_source)  ,mean)
AMG_reads_chayi_sum <- AMG_reads_chayi_sum[order(AMG_reads_chayi_sum$mapped_Normalized_Mapped_reads),]
AMG_reads_chayi_sum_a <- c("Carbohydrate metabolism","Glycan biosynthesis and metabolism","Biosynthesis of other secondary metabolites",
                           "Folding, sorting and degradation","Energy metabolism","Metabolism of cofactors and vitamins",
                           "Amino acid metabolism" ,"Metabolism of other amino acids")


############################## [Fig5-C] Painting #####################################
split.screen(c(length(unique(AMG_reads_chayi3$metabolism_source)), length(unique(AMG_reads_chayi3$mapped_Group)))) 
j=1
color <- c("#DCDDDD","#8A9CC4","#DA87B6","#8DC7B8","#C5B962")
names(color) <- unique(AMG_reads_chayi3$Bacteria_Host_G)
for (j in 1:length(unique(AMG_reads_chayi3$mapped_Group))){
  for (i in 1:length(unique(AMG_reads_chayi3$metabolism_source))){
    metabolism_source=AMG_reads_chayi_sum_a[i]
    mapped_Group=unique(AMG_reads_chayi3$mapped_Group)[j]
    group=paste(mapped_Group,metabolism_source,sep="_")
    AMG_reads_chayi4 <-AMG_reads_chayi3[AMG_reads_chayi3$Group2==group,]
    screen(6*i+j-6)
    if ( nrow(AMG_reads_chayi4)!=0  ){
      pie(as.numeric(AMG_reads_chayi4$mapped_Normalized_Mapped_reads),
          # labels =AMG_reads_chayi4$Bacteria_Host_G ,
          labels =NA ,
          # main = group,
          cex=0.3,cex.main=0.5,cex.axis=0.25,
          # radius =0.5,
          radius =log10(sum(AMG_reads_chayi4$mapped_Normalized_Mapped_reads)/sum(AMG_reads_chayi3$mapped_Normalized_Mapped_reads)+0.15)*1.2,
          col=color[AMG_reads_chayi4$Bacteria_Host_G],
          cex.lab=0.05 ,border = color[AMG_reads_chayi4$Bacteria_Host_G]
      )
    }
  }
  
}

close.screen(all = TRUE) 

ggsave("../fig5-C.pdf")
