##Script to Supplementary fig2-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
library()
library()

cran_packages=c("ggunchained","ggpubr","ggplot2")
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
col3_medians <- read.csv("../data/figS2C_1.csv")
ab <- read.csv("../data/figS2C_2.csv")

############################## Data processing ##############################
col3_medians$case <- col3_medians$diet
col3_medians$case[grep("Baseline",col3_medians$diet,ignore.case = T,)] <- "Baseline"
col3_medians$case[-grep("Baseline",col3_medians$diet,ignore.case = T,)] <- "Case"
col3_medians$Project <- sapply(strsplit(col3_medians$diet,split = "-"),"[",1)
col3_medians$Diet2 <- sapply(strsplit(col3_medians$diet,split = "-"),"[",2)
a <- unique(col3_medians[c("Project","Diet2")])
a <- a[duplicated(a$Project),]

col3_medians <- col3_medians[col3_medians$Project%in%a$Project,]
col3_medians <-merge(col3_medians,unique(ab[c("prj_name","type3")]),by.x="Project",by.y="prj_name")
a <- unique(col3_medians[c("Project","Diet2")])
a <- a[a$Diet2!="Baseline",]
a$diet_name <- paste0(a$Diet2,"_PRJ",1:nrow(a))
col3_medians <- merge(col3_medians,a[c("diet_name","Project")],by="Project")
col3_medians2 <- col3_medians
col3_medians2 <- col3_medians2[is.na(col3_medians2$dist)==F,]

a <- aggregate(list(dist=col3_medians2$dist),by=list(diet=col3_medians2$diet,
                                                     diet_name=col3_medians2$diet_name,
                                                     type3=col3_medians2$type3),median)
a <- a[-grep("Baseline",a$diet),]
a <- a[order(a$type3,a$dist,decreasing = T),]
col3_medians2$diet_name <- factor(col3_medians2$diet_name,levels=a$diet_name,labels=a$diet_name)
col3_medians2$dist[col3_medians2$dist<0.7] <- 0.7


############################## [FigS2-C] Painting #####################################
ggplot(col3_medians2,aes(x=diet_name,y=dist
                         ,fill=case))+
  geom_split_violin(trim = T,colour="white", scale = 'width')+
  facet_wrap(.~type3,scales = "free")+
  scale_fill_manual(values = c("#DEEBF7","#6BAED6"))+
  stat_compare_means(aes(group=case),method = "wilcox.test",label="p.signif",
                     show.legend = T,size=2)+
  # coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  stat_summary(fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               size=0.1,
               position = position_dodge(width = 0.3),
               alpha=0.5)

ggsave("../figS2-C.pdf") 
