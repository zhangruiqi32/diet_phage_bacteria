##Script to Figure 5-D and Supplementary fig8-B.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggalluvial","RColorBrewer","ggplot2")
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
meta <- read.table("~/cooperation/202409zhaofengxiang/data/sample_meta")
AMG_contig_contents <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_contig_contents.csv")
AMG_contents <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_contents.csv")
AMG_contig_length <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_contig_length.csv")
AMG_reads <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_reads.csv")
viruse_IMG_contig_anno <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_contig_anno.csv")

viruse_contig <- read.table("~/cooperation/202409zhaofengxiang/data/viruse_contig")
x <- read.csv("~/cooperation/202409zhaofengxiang/data/DABs_lefse.csv")


############################## Data processing ##############################
#AMG_reads<- merge(reads_number4s,meta,by="Sample_id")
#AMG_reads$Contig <- sapply(strsplit(AMG_reads$Contig,split ="-"),"[",1)
#AMG_reads$Contig2 <- substring(AMG_reads$Contig,1,30)
#AMG_contig_length <- read.table("/home/wanglab/database/diet/Fucoidan_5.4/forbwa/AMG_contig_length",sep="\t")
#colnames(AMG_contig_length) <- c("Contig","AMG_contig_length")

#AMG_reads<- merge(AMG_reads,AMG_contig_length,by.x="Contig",by.y="Contig2")
#AMG_reads$Normalized_Mapped_reads <- as.numeric(AMG_reads$reads_number)/(AMG_reads$AMG_contig_length*AMG_reads$num_seqs)*10000000000
#colnames(AMG_reads) <- paste("mapped",colnames(AMG_reads),sep="_")

#viruse_IMG_contig<- read.table("forzhushi/viruses_cd-hit2_blast_shai2",sep = "\t",header = F)
#viruse_IMG_contig_anno <- read.table("forzhushi/viruses_cd-hit2_blast_shai2_grep",header = F,sep = "\t")
#viruse_IMG_contig_anno$IMG_contig[viruse_IMG_contig_anno$V4!="whole"] <- paste(viruse_IMG_contig_anno$V1[viruse_IMG_contig_anno$V4!="whole"],viruse_IMG_contig_anno$V2[viruse_IMG_contig_anno$V4!="whole"],viruse_IMG_contig_anno$V3[viruse_IMG_contig_anno$V4!="whole"],viruse_IMG_contig_anno$V4[viruse_IMG_contig_anno$V4!="whole"],sep="|")
#viruse_IMG_contig_anno$IMG_contig[viruse_IMG_contig_anno$V4=="whole"] <- paste(viruse_IMG_contig_anno$V1[viruse_IMG_contig_anno$V4=="whole"],viruse_IMG_contig_anno$V2[viruse_IMG_contig_anno$V4=="whole"],viruse_IMG_contig_anno$V3[viruse_IMG_contig_anno$V4=="whole"],sep="|")
#viruse_IMG_contig_anno <- merge(viruse_IMG_contig[c(1,2)],viruse_IMG_contig_anno[c(5,15,17,20)],by.x="V2",by.y="IMG_contig",all.x=T)
#colnames(viruse_IMG_contig_anno) <- c("IMG_name","Contig_number","Source","Viruse_type","Bacteria_Host")

AMG_contig_length$Contig2 <- sapply(strsplit(AMG_contig_length$Contig,split ="-"),"[",1)
AMG_reads2_time  <- AMG_reads


viruse_IMG_contig_anno2<-viruse_IMG_contig_anno
viruse_IMG_contig_anno2$Bacteria_Host_G <-sapply(strsplit(viruse_IMG_contig_anno2$Bacteria_Host,split = "s__"),"[",1)
viruse_IMG_contig_anno2$Bacteria_Host_G <- sapply(strsplit(viruse_IMG_contig_anno2$Bacteria_Host_G,split = "g__"),"[",2)
viruse_IMG_contig_anno2$Bacteria_Host_G <- sapply(strsplit(viruse_IMG_contig_anno2$Bacteria_Host_G,split = ";"),"[",1)
viruse_IMG_contig_anno2$Contig_number2 <- substring(viruse_IMG_contig_anno2$Contig_number,1,30)
# AMG_reads2_time<- merge(AMG_reads2_time,viruse_IMG_contig_anno2,by.x="mapped_Contig2",by.y="Contig_number2")
viruse_contig <- read.table("/home/wanglab/database/diet/Fucoidan_5.4/viruse_contig")
AMG_reads2_time<- merge(AMG_reads2_time,viruse_contig,by.x="mapped_Contig2",by.y="Contig_number2")


AMG_reads3_time <- aggregate(list(mapped_Normalized_Mapped_reads=AMG_reads2_time$mapped_Normalized_Mapped_reads),
                             by=list(mapped_Group=AMG_reads2_time$mapped_Group,
                                     Contig_number=AMG_reads2_time$mapped_Contig,
                                     Bacteria_Host_G=AMG_reads2_time$Bacteria_Host_G,
                                     Bacteria_Host=AMG_reads2_time$Bacteria_Host),median)
# AMG_reads3_time <-AMG_reads2
AMG_reads3_time$Diet <- substring(AMG_reads3_time$mapped_Group,1,3)
AMG_reads3_time$Time<-sapply(strsplit(AMG_reads3_time$mapped_Group,split = "_"),"[",2)
AMG_reads3_time$Contig_number2 <- substring(AMG_reads3_time$Contig_number,1,30)

AMG_reads3_time <-merge(AMG_reads3_time,viruse_contig[c("Contig_number2","Phacts_result")],by.x="Contig_number2",by.y="Contig_number2",all.x=T)
# AMG_reads3_time$Bacteria_Host_G
# rongyuan_reads_forzhexian2 <- aggregate(rongyuan_reads_forzhexian$mapped_Normalized_Mapped_reads,by=list(mapped_Group=rongyuan_reads_forzhexian$mapped_Group,Contig_number=rongyuan_reads_forzhexian$mapped_Contig2),median )
AMG_reads3_time$sampleid_abundance <- paste(AMG_reads3_time$mapped_Group,AMG_reads3_time$Contig_number2,sep="_") 
# rongyuan_reads_forzhexian2$sampleid_abundance<- paste(rongyuan_reads_forzhexian2$mapped_Group,rongyuan_reads_forzhexian2$Contig_number,sep="_") 
# AMG_reads3_time <- merge(AMG_reads3_time,rongyuan_reads_forzhexian2[c("sampleid_abundance","x")],by="sampleid_abundance")



AMG_contig_contents$scaffold2 <- substring(AMG_contig_contents$scaffold, 1, last = 30)
AMG_contig_contents<-merge(AMG_contig_contents,meta,by="Sample_id")
AMG_contig_contents$Phacts_result[grep("Temperate|Lytic",AMG_contig_contents$Probability)]<- AMG_contig_contents$Probability[grep("Temperate|Lytic",AMG_contig_contents$Probability)]

AMG_contents2 <- aggregate(AMG_contents$Present.AMG.KOs,by=list(AMG_contents$Metabolism),function(x){paste(unique(x),collapse = ",")})
AMG_contents3 <- data.frame(matrix(ncol = 2,nrow = 0))
for (i in 1:nrow(AMG_contents2)){
  kos <- unlist(strsplit(AMG_contents2[i,2],split = ","))
  AMG_contents4 <-data.frame(metabolism=AMG_contents2[i,1],AMG.KO=kos)
  AMG_contents3 <- rbind(AMG_contents3,AMG_contents4)
}
AMG_contents3 <- unique(AMG_contents3)
AMG_contents_metabolism <-merge(AMG_contents3,AMG_contig_contents[,colnames(AMG_contig_contents)%in%grep("Sample_id",colnames(AMG_contig_contents),value = T)[-1]==F],by="AMG.KO")

AMG_contents2 <- aggregate(AMG_contents$Present.AMG.KOs,by=list(AMG_contents$Pathway),function(x){paste(unique(x),collapse = ",")})
AMG_contents3 <- data.frame(matrix(ncol = 2,nrow = 0))
for (i in 1:nrow(AMG_contents2)){
  kos <- unlist(strsplit(AMG_contents2[i,2],split = ","))
  AMG_contents4 <-data.frame(metabolism=AMG_contents2[i,1],AMG.KO=kos)
  AMG_contents3 <- rbind(AMG_contents3,AMG_contents4)
}
AMG_contents3 <- unique(AMG_contents3)
AMG_contents_Pathway <-merge(AMG_contents3,AMG_contig_contents[,colnames(AMG_contig_contents)%in%grep("Sample_id",colnames(AMG_contig_contents),value = T)[-1]==F],by="AMG.KO")


AMG_contents_metabolism2 <- AMG_contents_metabolism
colnames(AMG_contents_metabolism2) <- paste(colnames(AMG_contents_metabolism2),"source",sep="_")
AMG_contents_metabolism2<- unique(AMG_contents_metabolism2[c("protein_source","metabolism_source")])
AMG_reads3_time <- merge(AMG_reads3_time,AMG_contents_metabolism2,by.x="Contig_number",by.y="protein_source")

AMG_contents_metabolism2 <- AMG_contents_Pathway
colnames(AMG_contents_metabolism2) <- paste(colnames(AMG_contents_metabolism2),"source",sep="_")
AMG_contents_metabolism2<- unique(AMG_contents_metabolism2[c("protein_source","metabolism_source")])
colnames(AMG_contents_metabolism2) <- c("protein_source","pathway_source")
AMG_reads3_time <- merge(AMG_reads3_time,AMG_contents_metabolism2,by.x="Contig_number",by.y="protein_source")
length(table(AMG_reads3_time$pathway_source))
AMG_reads3_time$chayi <- "Normal bacteria"


chayijun <- sapply(strsplit(x$Spe,split = "g__"),"[",2)
chayijun <- paste(chayijun,collapse = "|")
AMG_reads3_time$chayi[grep(chayijun,AMG_reads3_time$Bacteria_Host)] <- "Differential bacteria"


############################## [FigS8-B] Painting #####################################
#diet-associated bacteria
a <- AMG_reads3_time[(AMG_reads3_time$metabolism_source=="Energy metabolism"|AMG_reads3_time$metabolism_source=="Amino acid metabolism")&
                       AMG_reads3_time$chayi=="Differential bacteria",]

a <- aggregate(list(mapped_Normalized_Mapped_reads=a$mapped_Normalized_Mapped_reads),
               by=list(Diet=a$Diet,Time=a$Time,
                       metabolism=a$metabolism_source),sum)

a$fill=paste(a$Diet,a$metabolism,sep="_")
a$Diet <- factor(a$Diet,levels = c('CON','HFD','FUC'),labels = c('CON','HFD','FUC'))
color<- c("#8BBCD6","#EE7C70","#ACD26A")
ggplot(a,aes(x=Time,y=mapped_Normalized_Mapped_reads,alluvium=Diet,stratum =Diet,fill=Diet,color=Diet ))+
  geom_flow(width = 0.3,#连线宽
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.2,#透明度
            color = 'white',#间隔颜色
            size = 0.5)+#间隔宽度
  geom_stratum(width = 0.28)+
  facet_wrap(~metabolism,scales = "free")+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+theme_test()


#all bacterial hosts
a <- AMG_reads3_time[(AMG_reads3_time$metabolism_source=="Carbohydrate metabolism"|AMG_reads3_time$metabolism_source=="Nucleotide metabolism"|
                        AMG_reads3_time$metabolism_source=="Lipid metabolism"),]
a <- aggregate(list(mapped_Normalized_Mapped_reads=a$mapped_Normalized_Mapped_reads),
               by=list(Diet=a$Diet,Time=a$Time,
                       metabolism=a$metabolism_source),sum)

a$fill=paste(a$Diet,a$metabolism,sep="_")
color<- c("#8BBCD6","#EE7C70","#ACD26A")
ggplot(a,aes(x=Time,y=mapped_Normalized_Mapped_reads,alluvium=Diet,stratum =Diet,fill=Diet,color=Diet ))+
  geom_flow(width = 0.3,#连线宽
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.2,#透明度
            color = 'white',#间隔颜色
            size = 0.5)+#间隔宽度
  geom_stratum(width = 0.28)+
  facet_wrap(~metabolism,scales = "free")+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+theme_test()


############################## [Fig5-D] Painting #####################################
metabolism <- "Amino acid metabolism"

a <- AMG_reads3_time[(AMG_reads3_time$metabolism_source==metabolism)&
                       AMG_reads3_time$chayi=="Differential bacteria",]

a$Contig_number2 <- substring(a$Contig_number,1,30)
a <- AMG_reads2_time[AMG_reads2_time$mapped_Contig%in%a$Contig_number,]
a <- merge(a,unique(AMG_reads3_time[c("Contig_number","metabolism_source" ,"pathway_source" )]),by.y="Contig_number",by.x="mapped_Contig")
a <- a[a$metabolism_source==metabolism,]
a<- aggregate(list(AMG_Abundance=a$mapped_Normalized_Mapped_reads),
              by=list(Mapped_GD=a$mapped_GD,
                      Mapped_Diet=a$mapped_Diet,
                      Bacteria_Host_G=a$Bacteria_Host_G,
                      metabolism_source=a$metabolism_source,
                      Sample_id=a$mapped_Sample_id) ,sum
)

a1<- a
a$AMG_Abundance<- log10(a$AMG_Abundance)

median <- aggregate(list(AMG_Abundance=a1$AMG_Abundance),by=list(Mapped_Diet=a1$Mapped_Diet,
                                                           Bacteria_Host_G=a1$Bacteria_Host_G,
                                                           Mapped_GD=a1$Mapped_GD) ,mean)
median$AMG_Abundance <- log10(median$AMG_Abundance)
sd <-  aggregate(list(max=a$AMG_Abundance),by=list(Mapped_Diet=a$Mapped_Diet,
                                                Bacteria_Host_G=a$Bacteria_Host_G,
                                                Mapped_GD=a$Mapped_GD) ,sd)

sd <-  aggregate(list(max=a$AMG_Abundance),by=list(Mapped_Diet=a$Mapped_Diet,
                                                Bacteria_Host_G=a$Bacteria_Host_G,
                                                Mapped_GD=a$Mapped_GD) ,function(x)
                                                  sd(x)/length(x) )

a <- data.frame(median,sd$max)
a$Mapped_Diet <- factor(a$Mapped_Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC") )
a$Bacteria_Host_G <- paste(a$Bacteria_Host_G,"Phage",sep="_")

ggplot(a,aes(x=Mapped_GD,y=AMG_Abundance,group=Bacteria_Host_G,color=Bacteria_Host_G))+
  geom_point(size=3,alpha=0.8,shape=19)+facet_wrap(.~Mapped_Diet)+
  geom_line(position = position_dodge(0.1),cex=1.3,alpha=0.8)+xlab("Time")+ylab("Log10(AMG Abundance+1)")+
  geom_errorbar(aes(x =Mapped_GD,ymin = AMG_Abundance-sd.max, ymax =AMG_Abundance+sd.max,width=.2))+
  scale_color_manual(values = c("#DCDDDD","#8DA0CB","#E78AC3")  )+theme_classic()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background = element_rect(
          # color="gray", linetype="solid"
          color="white",fill="white", linetype="solid"
        ),legend.position = "bottom",
        axis.title=element_text(size=11),legend.title = element_text(size=11)
  )
ggsave("../fig5-D.pdf")


