##Script to Supplementary fig8-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","ggsignif","reshape2","ggpubr","ggTimeSeries","RColorBrewer","ggplot2")
bioconductor_packages=c()

for(package in cran_packages){
    library(package, character.only=T,quietly = T, warn.conflicts = F)
}

if (!require("BiocManager", quietly = T, warn.conflicts = F)) install.packages("BiocManager")

for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

library(tidyverse)
library(reshape2)
############################## Data loading ##############################
HGT_reads <- read.csv("~/cooperation/202409zhaofengxiang/data/fig4BC_S8A_reads.csv")
HGT_contig_length <- read.csv("~/cooperation/202409zhaofengxiang/data/fig4BC_S8A_length.csv")

viruse_contig <- read.table("~/cooperation/202409zhaofengxiang/data/viruse_contig")
x <- read.csv("~/cooperation/202409zhaofengxiang/data/DABs_lefse.csv")
meta <- read.table("~/cooperation/202409zhaofengxiang/data/sample_meta")



reads<- read.table("~/cooperation/202409zhaofengxiang/data/figS9A_long_format_stats.tsv",header = T)
AMGs_length <- read.csv("~/cooperation/202409zhaofengxiang/data/figS9A_combined_sequences_lengths.txt",
                        header =F,sep = "\t",col.names = c("Contig","Length"))

############################## Data processing ##############################
HGT_contig_length$Contig2 <- sapply(strsplit(HGT_contig_length$Contig,split ="-"),"[",1)
HGT_reads<- merge(HGT_reads,HGT_contig_length,by.x="Contig",by.y="Contig2")
HGT_reads$Normalized_Mapped_reads <- as.numeric(HGT_reads$reads_number)/(as.numeric(HGT_reads$HGT_contig_length)*as.numeric(HGT_reads$num_seqs))*10000000000
HGT_reads2 <-spread(HGT_reads[colnames(HGT_reads)%in%c("num_seqs","reads_number","GD","Diet","Group")==F],key=Sample_id,value =Normalized_Mapped_reads )
HGT_reads2[is.na(HGT_reads2)] <- 0
HGT_reads2 <-gather(HGT_reads2,key = "Sample_id",value = "Normalized_Mapped_reads",L1EGH020932:L1HGH010566)
HGT_reads2 <-merge(HGT_reads2,meta,by="Sample_id")

files <- paste0("/wang/students/zhaofengxiang/Fucoidan_5.4/",unique(meta$Sample_id),"/")
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
HGT_contig_contents <- HGT_contig_contents[grep("integrase|recombinase|transposase|holins|lysins|spanins|repressor",HGT_contig_contents$V5,ignore.case = T),]

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
HGT_contig_contents2 <- HGT_contig_contents2[is.null(HGT_contig_contents2$type)==F,]

chayijun <- sapply(strsplit(x$Spe,split = "g__"),"[",2)
chayijun <- paste(chayijun,collapse = "|")
viruse_IMG_contig_anno_chayijun <- viruse_contig[grep(chayijun,viruse_contig$Bacteria_Host),]

viruse_IMG_contig_anno_chayijun$Contig_name2 <- viruse_IMG_contig_anno_chayijun$Contig_number
viruse_IMG_contig_anno_chayijun$Contig_name2 <- substring(viruse_IMG_contig_anno_chayijun$Contig_name2, 1, last = 30)
HGT_viruse_IMG_contig_anno_chayijun <- data.frame(viruse_IMG_contig_anno_chayijun[viruse_IMG_contig_anno_chayijun$Contig_name2%in%HGT_contig_contents2$scaffold2,])

# a<- HGT_viruse_IMG_contig_anno_chayijun
# write.table(HGT_viruse_IMG_contig_anno_chayijun$Contig_number,"/wang/students/zhaofengxiang/Fucoidan_5.4/depth/contigs_name_integrase",
#             row.names = F,col.names = F,quote = F)

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

##添加HGT相关基因
HGT_list <- lapply(1:length(HGT_contigs), function(i) {
  df <- read.csv(HGT_contigs[i], sep = "\t", header = FALSE)
  if (nrow(df) > 0) df$Sample_id <- files[i]
  df
})

HGT_contig_contents2 <- do.call(rbind, HGT_list)
HGT_contig_contents2 <- subset(HGT_contig_contents2, 
                               grepl(">", V1) & grepl("integrase|recombinase|transposase|holin|lysin|spanins|repressor", V5, ignore.case = TRUE))
HGT_contig_contents2$scaffold2 <- substring(HGT_contig_contents2$V1, 2, 31)
HGT_contig_contents2$type <- sapply(HGT_contig_contents2$V5, function(desc) {
  desc_lower <- tolower(desc)
  patterns <- c("integrase", "recombinase", "transposase", 
                "holin", "lysin", "spanins", "repressor")
  matches <- sapply(patterns, function(p) grepl(p, desc_lower))
  if (any(matches)) patterns[which(matches)[1]] else NA
})
HGT_contig_contents2$gene <- gsub(">","",HGT_contig_contents2$V1)
HGT_contig_contents2$gene <- gsub(" ","_",HGT_contig_contents2$gene)
HGT_contig_contents3<- HGT_contig_contents2
#write.table(HGT_contig_contents2$gene,"/wang/students/zhaofengxiang/Fucoidan_5.4/12.8_new/HGT_genes/lifestyle_contigs_id.txt",col.names = F,row.names = F,quote = F)

reads <- merge(reads,meta,by.x = "Sample",by.y="Sample_id")
reads <- merge(reads,AMGs_length,by.x = "Contig",by.y="Contig")
reads$Normalized_Mapped_reads <- reads$Reads/(as.numeric(reads$num_seqs)*as.numeric(reads$Length))*10^9
reads2 <- reads[-grep("bin",reads$Contig),]
reads2 <- merge(reads2,HGT_contig_contents3,by.x="Contig",by.y = "gene")
# reads2 <- reads2[reads2$Contig%in%HGT_viruse_IMG_contig_anno2$Contig,]
unique(HGT_viruse_IMG_contig_anno2$Contig[HGT_viruse_IMG_contig_anno2$Contig%in%reads2$Contig==F])

plot_data <- reads2[reads2$Normalized_Mapped_reads!=0,] %>%
  group_by(type, Sample, Diet, GD) %>%
  summarise(expr = sum(Normalized_Mapped_reads), .groups = "drop") %>%
  group_by(type, Diet, GD) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE), #median
    sd_expr = sd(expr, na.rm = TRUE),
    se_expr = sd_expr / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(type, Diet, GD)


plot_data_log2 <- plot_data %>%
   mutate( mean_expr_log2= log2(mean_expr + 1),
           se_expr_log2= sd_expr / (sqrt(n_samples) * (mean_expr + 1) * log(2)))
  #mutate( mean_expr_log2= mean_expr,
  #      se_expr_log2= sd_expr / (sqrt(n_samples) * (mean_expr + 1) ))
# plot_data_log2 <- plot_data_log2 %>%filter(Contig %in% result$Contig   )
plot_data_log2$Diet <- factor(plot_data_log2$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC") )


############################# [FigS8-C] Painting #####################################
# DIET-GENE-log2分面折线图
p1<- ggplot(plot_data_log2, aes(x = GD, y = mean_expr_log2, color = Diet, group = Diet)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_expr_log2 - se_expr_log2, ymax = mean_expr_log2 + se_expr_log2), 
                width = 0.1, size = 0.5) +
  facet_grid(rows = vars(type), cols = vars(Diet), scales = "free_y") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  # 坐标轴标签
  labs(
    x = "Time (Day)",
    y = "log2(Normalized Expression)",
    color = "Diet",
    title = "HGT-related genes"
  ) +
  # 颜色设置
  scale_color_manual(values = c(CON="#8BBCD6","HFD" = "#EE7C70", "FUC" = "#ACD26A"))
p1

# GENE分面折线图
plot_data_figS9A <- plot_data %>%
  filter(type %in% c("holin", "lysin"))
plot_data_figS9A$Diet <- factor(plot_data_figS9A$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC") )

ggplot(plot_data_figS9A, aes(x = GD, y = mean_expr, color = Diet, group = Diet)) +
  geom_line(position = position_dodge(0.1),cex=1.3,alpha=0.9) +
  geom_point(size=3,alpha=0.9,shape=19) +
  geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr), 
                width = 0.1, size = 0.3, alpha=0.5) +
  facet_grid(rows = vars(type), 
             #cols = vars(Diet), 
             scales = "free_y") +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  # 坐标轴标签
  labs(
    x = "Time (Day)",
    y = "Normalized Expression",
    color = "Diet",
    title = "lysis-related genes"
  ) +
  # 颜色设置
  scale_color_manual(values = c(CON="#DCDDDD","HFD" = "#EE7C70", "FUC" = "#ACD26A"))

