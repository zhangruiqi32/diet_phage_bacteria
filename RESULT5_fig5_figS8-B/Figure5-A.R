##Script to Figure 5-A.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("dplyr","tidyr","ggraph","igraph","ggnewscale","RColorBrewer","ggplot2")
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
meta <- read.table("../data/sample_meta")
AMG_contig_contents <- read.csv("../data/fig5_AMG_contig_contents.csv")
AMG_contents <- read.csv("../data/fig5_AMG_contents.csv")
AMG_contig_length <- read.csv("../data/fig5_contig_length.csv")
AMG_reads <- read.csv("../data/fig5_AMG_reads.csv")
viruse_IMG_contig_anno <- read.csv("../data/fig5_contig_anno.csv")

viruse_contig <- read.table("../data/viruse_contig")
x <- read.csv("../data/DABs_lefse.csv")

############################## Data processing ##############################
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
#viruse_IMG_contig_anno<-merge(viruse_IMG_contig[c(1,2)],viruse_IMG_contig_anno[c(5,15,17,20)],by.x="V2",by.y="IMG_contig",all.x=T)
#colnames(viruse_IMG_contig_anno) <- c("IMG_name","Contig_number","Source","Viruse_type","Bacteria_Host")

AMG_contig_length$Contig2 <- sapply(strsplit(AMG_contig_length$Contig,split ="-"),"[",1)
AMG_reads2  <- AMG_reads[grep("D0|D10",AMG_reads$mapped_GD),]

viruse_IMG_contig_anno2 <- viruse_IMG_contig_anno
viruse_IMG_contig_anno2$Contig_number2 <- substring(viruse_IMG_contig_anno2$Contig_number,1,30)
viruse_IMG_contig_anno2$Bacteria_Host_G <- sapply(strsplit(viruse_IMG_contig_anno2$Bacteria_Host,split = ";g__" ) ,"[",2)
viruse_IMG_contig_anno2$Bacteria_Host_G <- sapply(strsplit(viruse_IMG_contig_anno2$Bacteria_Host_G,split = ";" ) ,"[",1)
AMG_reads2<- merge(AMG_reads2,viruse_IMG_contig_anno2,by.x="mapped_Contig2",by.y="Contig_number2")
AMG_reads2$fill <- paste(AMG_reads2$Bacteria_Host_G,AMG_reads2$mapped_GD,sep="_")

AMG_reads3 <- aggregate(list(mapped_Normalized_Mapped_reads=AMG_reads2$mapped_Normalized_Mapped_reads),by=list(mapped_Group=AMG_reads2$mapped_Group,Contig_number=AMG_reads2$mapped_Contig),median)
AMG_reads3$Diet <- sapply(strsplit(AMG_reads3$mapped_Group,split = "_"),"[",1)
AMG_reads3$Time<- sapply(strsplit(AMG_reads3$mapped_Group,split = "_"),"[",2)
AMG_reads3$Contig_number2 <- substring(AMG_reads3$Contig_number,1,30)
AMG_reads3 <-merge(AMG_reads3,viruse_IMG_contig_anno2[c("Contig_number2","Bacteria_Host_G","Bacteria_Host")],by.x="Contig_number2",by.y="Contig_number2")
AMG_reads3 <- merge(AMG_reads3,viruse_contig[c("Contig_number2","Phacts_result")],by="Contig_number2")

AMG_contents_metabolism2 <- AMG_contents_metabolism
colnames(AMG_contents_metabolism2) <- paste(colnames(AMG_contents_metabolism2),"source",sep="_")
AMG_contents_metabolism2<- unique(AMG_contents_metabolism2[c("protein_source","metabolism_source")])
AMG_reads3 <- merge(AMG_reads3,AMG_contents_metabolism2,by.x="Contig_number",by.y="protein_source")

AMG_contents_metabolism2 <- AMG_contents_Pathway
colnames(AMG_contents_metabolism2) <- paste(colnames(AMG_contents_metabolism2),"source",sep="_")
AMG_contents_metabolism2<- unique(AMG_contents_metabolism2[c("protein_source","metabolism_source")])
colnames(AMG_contents_metabolism2) <- c("protein_source","pathway_source")
AMG_reads3 <- merge(AMG_reads3,AMG_contents_metabolism2,by.x="Contig_number",by.y="protein_source")
length(table(AMG_reads3$pathway_source))


chayijun <- sapply(strsplit(x$Spe,split = "g__"),"[",2)
chayijun <- paste(chayijun,collapse = "|")
AMG_reads3$chayi <- "Normal bacteria"
AMG_reads3$chayi[grep(chayijun,AMG_reads3$Bacteria_Host)] <- "Differential bacteria"

AMG_reads4 <-spread(AMG_reads3[colnames(AMG_reads3)%in%c("mapped_Group")==F], key= 'Time', value =  'mapped_Normalized_Mapped_reads')
AMG_reads4 <- aggregate(list(D0=AMG_reads4$D0,D10=AMG_reads4$D10),
                        by=list(Contig_number=AMG_reads4$Contig_number,metabolism=AMG_reads4$metabolism_source,pathway=AMG_reads4$pathway_source,chayi=AMG_reads4$chayi),sum)
AMG_reads4$D10[is.na(AMG_reads4$D10)] <- 0
AMG_reads4$D0[is.na(AMG_reads4$D0)] <- 0
AMG_reads4$sum <- AMG_reads4$D0+AMG_reads4$D10
AMG_reads4$sum[AMG_reads4$chayi=="Differential bacteria"] <- AMG_reads4$sum[AMG_reads4$chayi=="Differential bacteria"]/sum(AMG_reads4$sum[AMG_reads4$chayi=="Differential bacteria"])
AMG_reads4$sum[AMG_reads4$chayi!="Differential bacteria"] <- AMG_reads4$sum[AMG_reads4$chayi!="Differential bacteria"]/sum(AMG_reads4$sum[AMG_reads4$chayi!="Differential bacteria"])
AMG_reads4_diff <- AMG_reads4[AMG_reads4$chayi=="Differential bacteria",]
AMG_reads4_diff_metab <- aggregate(list(sum=AMG_reads4_diff$sum),by=list(metabolism=AMG_reads4_diff$metabolism),sum)
AMG_reads4_diff_metab <- AMG_reads4_diff_metab[order(AMG_reads4_diff_metab$sum),]
AMG_reads4_diff_metab$cm <- paste("Differential bacteria",AMG_reads4_diff_metab$metabolism,sep="_")
AMG_reads4_diff$metabolism <- factor(AMG_reads4_diff$metabolism,levels = AMG_reads4_diff_metab$metabolism,labels = AMG_reads4_diff_metab$metabolism)
AMG_reads4_diff_pathway <- aggregate(list(sum=AMG_reads4_diff$sum),by=list(metabolism=AMG_reads4_diff$metabolism,pathway=AMG_reads4_diff$pathway),sum)
AMG_reads4_diff_pathway$metabolism <- factor(AMG_reads4_diff_pathway$metabolism,levels = AMG_reads4_diff_metab$metabolism,labels = AMG_reads4_diff_metab$metabolism)
AMG_reads4_diff_pathway$pathway[AMG_reads4_diff_pathway$sum<0.02] <- ""
AMG_reads4_diff_pathway <- AMG_reads4_diff_pathway[order(AMG_reads4_diff_pathway$metabolism,AMG_reads4_diff_pathway$sum),]
AMG_reads4_diff_pathway$mp <- paste(AMG_reads4_diff_pathway$metabolism,AMG_reads4_diff_pathway$pathway,sep="_")
AMG_reads4_diff_pathway$cmp <- paste("Differential bacteria",AMG_reads4_diff_pathway$mp,sep="_")
AMG_reads4_diff$mp <- paste(AMG_reads4_diff$metabolism,AMG_reads4_diff$pathway,sep="_")
AMG_reads4_diff$mp <- factor(AMG_reads4_diff$mp,levels = AMG_reads4_diff_pathway$mp,labels = AMG_reads4_diff_pathway$mp)
AMG_reads4_diff <- AMG_reads4_diff[order(AMG_reads4_diff$metabolism,AMG_reads4_diff$mp),]

AMG_reads4_normal <- AMG_reads4[AMG_reads4$chayi!="Differential bacteria",]
AMG_reads4_normal_metab <- aggregate(list(sum=AMG_reads4_normal$sum),by=list(metabolism=AMG_reads4_normal$metabolism),sum)
AMG_reads4_normal_metab <- AMG_reads4_normal_metab[order(AMG_reads4_normal_metab$sum,decreasing = T),]
AMG_reads4_normal_metab$cm <- paste("Normal bacteria",AMG_reads4_normal_metab$metabolism,sep="_")
AMG_reads4_normal$metabolism <- factor(AMG_reads4_normal$metabolism,levels = AMG_reads4_normal_metab$metabolism,labels = AMG_reads4_normal_metab$metabolism)
AMG_reads4_normal_pathway <- aggregate(list(sum=AMG_reads4_normal$sum),by=list(metabolism=AMG_reads4_normal$metabolism,pathway=AMG_reads4_normal$pathway),sum)
# AMG_reads4_normal_pathway$pathway[AMG_reads4_normal_pathway$sum<0.02] <- "Others"
AMG_reads4_normal_pathway$pathway[AMG_reads4_normal_pathway$sum<0.02] <- ""
AMG_reads4_normal_pathway <-aggregate(list(sum=AMG_reads4_normal_pathway$sum),by=list(metabolism=AMG_reads4_normal_pathway$metabolism,
                                                                                      pathway=AMG_reads4_normal_pathway$pathway      ), sum ) 
AMG_reads4_normal_pathway$metabolism <- factor(AMG_reads4_normal_pathway$metabolism,levels = AMG_reads4_normal_metab$metabolism,labels = AMG_reads4_normal_metab$metabolism)
AMG_reads4_normal_pathway <- AMG_reads4_normal_pathway[order(AMG_reads4_normal_pathway$metabolism,-AMG_reads4_normal_pathway$sum),]

AMG_reads4_normal_pathway$mp <- paste(AMG_reads4_normal_pathway$metabolism,AMG_reads4_normal_pathway$pathway,sep="_")
AMG_reads4_normal_pathway$cmp <- paste("Normal bacteria",AMG_reads4_normal_pathway$mp,sep="_")
AMG_reads4_normal$mp <- paste(AMG_reads4_normal$metabolism,AMG_reads4_normal$pathway,sep="_")
AMG_reads4_normal$mp <- factor(AMG_reads4_normal$mp,levels = AMG_reads4_normal_pathway$mp,labels = AMG_reads4_normal_pathway$mp)
AMG_reads4_normal <- AMG_reads4_normal[order(AMG_reads4_normal$metabolism,AMG_reads4_normal$mp),]


AMG_reads4_diff_metab$chayi <- "Differential bacteria"
AMG_reads4_normal_metab$chayi <- "Normal bacteria"
AMG_reads4_diff_pathway$cm <- paste("Differential bacteria",AMG_reads4_diff_pathway$metabolism,sep="_")
AMG_reads4_normal_pathway$cm <- paste("Normal bacteria",AMG_reads4_normal_pathway$metabolism,sep="_")
edges<-data.frame(rbind(
  cbind(rep('origin',2), c("Normal bacteria","Differential bacteria") ),  
  as.matrix(unique(AMG_reads4_normal_metab[c("chayi","cm")])),
  as.matrix(unique(AMG_reads4_normal_pathway[c("cm","cmp")])), 
  as.matrix(unique(AMG_reads4_diff_metab[c("chayi","cm")])),
  as.matrix(unique(AMG_reads4_diff_pathway[c("cm","cmp")]))   ) )

colnames(edges)<-c('from','to')

vertices0<-data.frame(name=unique(c(as.character(edges$from), as.character(edges$to))))
df_leaf<-data.frame(rbind(
  cbind(unique(as.character(AMG_reads4$chayi)),rep(0.1,2), unique(as.character(AMG_reads4$chayi)) ),
  as.matrix(unique(AMG_reads4_diff_metab[c("cm","sum","metabolism")])),
  as.matrix(unique(AMG_reads4_diff_pathway[c("cmp","sum","pathway")])),
  as.matrix(unique(AMG_reads4_normal_metab[c("cm","sum","metabolism")])),
  as.matrix(unique(AMG_reads4_normal_pathway[c("cmp","sum","pathway")]))  ) )
colnames(df_leaf) <- c('cmp','Value','label')
df_leaf$angle<-90-(1:nrow(df_leaf))/nrow(df_leaf)*360
df_leaf$angle<-ifelse(df_leaf$angle< -90, df_leaf$angle+180, df_leaf$angle)
vertices<-left_join(vertices0,df_leaf,by=c('name'='cmp'))
vertices$name <- factor(vertices$name,levels = vertices$name,labels = vertices$name)
vertices<-vertices[order(vertices$name),]
vertices$name <- as.character(vertices$name)
vertices$name2 <- sapply(strsplit(vertices$name,split = "_"),"[",2)


############################## [Fig5-A] Painting #####################################
color_metab <-c("#BCB7D6","#89CCC0","#FAE4CE","#ADD26A",
                "#EE7C6F","#7FABCB","#F8F5B3","#F6CBDE",
                "#DCDDDD","#DCDDDD","#DCDDDD","#DCDDDD")
names(color_metab)<- unique(c(AMG_reads4_diff_metab$metabolism,AMG_reads4_normal_metab$metabolism))

graph <- graph_from_data_frame(edges, vertices = vertices)
ggraph(graph, layout ='partition', circular = TRUE,weight=as.numeric(Value)) + 
  geom_node_arc_bar(aes(filter =(depth==1 ),
                        colour=label,r0=0.3,r=0.6),
                    fill=NA,size=0.5)+
  scale_colour_manual(values = c("Normal bacteria" = "#1E2C5B", 
                                 "Differential bacteria" = "#A92424"))+
  new_scale_colour()+
  geom_node_arc_bar(aes(filter =(depth==2 ),
                        fill=label,r0=0.6,r=0.9),
                    colour="#727171",size=0.3,alpha=0.5)+
  scale_fill_manual(values = color_metab)+ 
  new_scale_fill()+
  geom_node_arc_bar(aes(filter =(depth==3 ),
                        fill=name2,r0=0.9,r=1.5),
                    colour="#727171",size=0.3,alpha=0.5)+
  scale_fill_manual(values = color_metab)+
  geom_node_text(aes(filter =(depth==3),
                     label=label,angle=0),
                 size=2,colour="black",
                 hjust ="inward") +
  coord_fixed()+
  guides(size="none",fill="none")+
  theme_void()
ggsave("../fig5-A.pdf")
