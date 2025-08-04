##Script to Figure 1-F-3 and Supplementary fig4.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","RColorBrewer","ggplot2")
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
dfs_crass2 <- read.csv("../data/fig1F-3_1.csv")
crass_phage2 <- read.csv("../data/fig1F-3_2.csv")

prj2 <- read.csv("../data/fig1F_prj.csv")
diet_controls2 <- read.csv("../data/fig1F_diet2.csv")

crass_phage3 <- read.table("../data/fig1_crass_phage3")

############################## [FigS4] Processing and Painting #####################################
colnames(dfs_crass2) <- c("Bact","Viruse","p_value","R2","PRJ-DIET")    
dfs_crass2 <- dfs_crass2[-grep("_nutrition",dfs_crass2$`PRJ-DIET`,fixed=T),]
dfs_crass2$phylum <- sapply(strsplit(dfs_crass2$Bact,split="g__",fixed = T ),"[",2)
dfs_crass2$phylum[grep("crass",dfs_crass2$Viruse,ignore.case = T)] <- "crAssphage"
dfs_crass2$phylum <- sapply(strsplit(dfs_crass2$phylum,split="|",fixed = T ),"[",1)
dfs_xiangguanxinngs <- dfs_crass2

######控制去p值
dfs_xiangguanxinngs <- dfs_xiangguanxinngs[dfs_xiangguanxinngs$p_value<0.05,]

dfs_xiangguanxinngs$Diet <- sapply(strsplit(dfs_xiangguanxinngs$`PRJ-DIET`,split="-"),"[",2)
dfs_xiangguanxinngs$Diet <- gsub("/","",dfs_xiangguanxinngs$Diet,fixed = T)

dfs_xiangguanxinngs$Project <- sapply(strsplit(dfs_xiangguanxinngs$`PRJ-DIET`,split="-"),"[",1)
dfs_xiangguanxinngs <- dfs_xiangguanxinngs[dfs_xiangguanxinngs$Project%in%unique(crass_phage2$Project),]


i=1
dfs_xiangguanxinngs_CTs <- data.frame(matrix(ncol=1,nrow = 0))
for (i in 1:nrow(diet_controls2) ){
  diet=diet_controls2$case[i]
  diet_control <- unlist(strsplit(diet_controls2$control[i],split = "/"))
  
  prjs <- unique(prj2$BioProject[prj2$Diet1==diet])
  j=1
  for (j in 1:length(prjs)){
    dfs_xiangguanxinngs_CT <- dfs_xiangguanxinngs[dfs_xiangguanxinngs$Project%in%prjs[j],]
    dfs_xiangguanxinngs_CT<- dfs_xiangguanxinngs_CT[dfs_xiangguanxinngs_CT$Diet%in%c(diet_control,diet),]
    
    if (nrow(dfs_xiangguanxinngs_CT)!=0){
      dfs_xiangguanxinngs_treat <- dfs_xiangguanxinngs_CT[dfs_xiangguanxinngs_CT$Diet%in%diet,]
      if(diet_control[1]%in%dfs_xiangguanxinngs_CT$Diet ){
        dfs_xiangguanxinngs_control <- dfs_xiangguanxinngs_CT[dfs_xiangguanxinngs_CT$Diet%in%diet_control[1],]
      }else{
        dfs_xiangguanxinngs_control <- dfs_xiangguanxinngs_CT[dfs_xiangguanxinngs_CT$Diet%in%diet_control[2],]
      }
      
      dfs_xiangguanxinngs_CTs <- rbind(dfs_xiangguanxinngs_CTs,
                                       data.frame(dfs_xiangguanxinngs_treat,type="case",prj_case_diet=paste(prjs[j],diet,sep="-")   ),
                                       data.frame(dfs_xiangguanxinngs_control,type="control",prj_case_diet=paste(prjs[j],diet,sep="-")   )
      )
    } 
  }
}
dfs_xiangguanxinngs<- dfs_xiangguanxinngs_CTs

dfs_xiangguanxinngs$Diet <- sapply(strsplit(dfs_xiangguanxinngs$prj_case_diet,split="-"),"[",2)
dfs_xiangguanxinngs$Diet2 <- dfs_xiangguanxinngs$Diet
dfs_xiangguanxinngs$Diet2 <- gsub("_"," ",dfs_xiangguanxinngs$Diet2)
dfs_xiangguanxinngs$Diet2 <- gsub("diet","",dfs_xiangguanxinngs$Diet2)
dfs_xiangguanxinngs$Diet <- factor(dfs_xiangguanxinngs$Diet,levels = unique(dfs_xiangguanxinngs$Diet),
                                   labels = unique(dfs_xiangguanxinngs$Diet2))

crass_phage2<- dfs_xiangguanxinngs
crass_phage2$x <- paste(crass_phage2$phylum,crass_phage2$Diet,sep="_")

unique(crass_phage3$Diet)
crass_phage2<- crass_phage2[crass_phage2$phylum%in%unique(crass_phage3$phylum),]
crass_phage2$phylum <- factor(crass_phage2$phylum,levels =unique(crass_phage3$phylum) ,labels = unique(crass_phage3$phylum))
crass_phage2<- crass_phage2[crass_phage2$Diet%in%unique(crass_phage3$Diet),]
crass_phage2$Diet <- factor(crass_phage2$Diet,levels =unique(crass_phage3$Diet) ,labels = unique(crass_phage3$Diet))
crass_phage2<- crass_phage2[crass_phage2$Project%in%unique(crass_phage3$Project),]

crass_phage2$R2 <- abs(crass_phage2$R2)
crass_phage3 <- crass_phage2
crass_phage3 <- crass_phage3[order(crass_phage3$phylum ,crass_phage3$Diet,crass_phage3$R2),]
crass_phage3$x2 <- paste(crass_phage3$phylum ,crass_phage3$Diet,crass_phage3$type,sep="_")
crass_phage3$x2 <- factor(crass_phage3$x2,levels = unique(crass_phage3$x2),labels = unique(crass_phage3$x2))
crass_phage3 <- crass_phage3[order(crass_phage3$x2,crass_phage3$R2),]

a<- data.frame(unlist(aggregate(1:nrow(crass_phage3),by=list(x2=crass_phage3$x2),function(x) 1:length(x) )[2] ) )
crass_phage3$number <- unlist(a[1])
crass_phage3$number <- factor(crass_phage3$number,levels = crass_phage3$number,labels = crass_phage3$number)

##### Painting
c <- c(  rep("top",length(unique(crass_phage3$Diet))),  rep("right",length(unique(crass_phage3$Diet)))  )
ggplot(crass_phage3,aes(x=number,y=R2,color=type))+
  geom_point(shape=16,size=0.5,alpha=0.8)+
  facet_wrap(phylum~Diet,scales = "free_x",shrink = T,as.table = T,drop =F,ncol = length(unique(crass_phage3$Diet)) ) +
  #facet_grid(phylum~Diet,scales = "free_x",shrink = T,as.table = T,drop =F) +
  scale_color_manual(values = c(brewer.pal(10,"Paired")[5],"gray"))+
  # scale_color_manual(values = a)+
  # coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "#727171", size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "#727171", linewidth = 0.25,lineend = "square"),
        axis.ticks.length.y =  unit(0.25,"mm"),
        strip.text = element_blank(),
        #strip.text = element_text(size = 5), strip.background = element_blank(),
        panel.grid=element_blank() ,
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(color = "#727171", linewidth = 0.5, fill = NA),
        panel.spacing.x = unit(1.2, "mm"),
        panel.spacing.y = unit(1.8, "mm"))
ggsave("../figS4.pdf") 


############################## [Fig1-F-3] Processing and Painting #####################################
#### by list按crass_phage3中的菌属
##以噬菌体（s__）——phylum
R_median <- aggregate(list(median=abs(crass_phage3$R2)),
                      by=list(Diet=crass_phage3$Diet,
                              phylum=crass_phage3$phylum,
                              type=crass_phage3$type),
                      median) #数据框crass_phage3中，按照phylum、Diet和type这三个变量进行分组后，每组中R2列绝对值的中位数
R_median_wide <- tidyr::spread(R_median, key=type, value=median ) #存在缺失值
R_pvalue_rawdata <- dplyr::left_join(crass_phage3[,4:12],R_median_wide, by=c("Diet","phylum")) #crass_phage3要包含type和x那列
R_pvalue_rawdata <- na.omit(R_pvalue_rawdata) 
R_pvalue_rawdata$type <- as.factor(R_pvalue_rawdata$type)
R_pvalue_rawdata$type <- droplevels(R_pvalue_rawdata$type)
R_pvalue_rawdata$x <- as.factor(R_pvalue_rawdata$x)
R_pvalue_rawdata$x <- droplevels(R_pvalue_rawdata$x)
R_pvalue <- aggregate(list(pvalue=1:nrow(R_pvalue_rawdata)),
                      by=list(phylum_Diet=R_pvalue_rawdata$x), #有的组只有case/control 不能wilcox 需要删掉
                      function(x)wilcox.test(R2~type,R_pvalue_rawdata[x,])[["p.value"]])#类似箱线图 wilcox检验显著性 但判读上升/下降 是比两组的中位值
R_median_wide$aftertreat <- ifelse(R_median_wide$case>R_median_wide$control, "up", "down") #菌与噬菌体相关性：饮食干预后上升
R_median_wide$phylum_Diet <- paste(R_median_wide$phylum, R_median_wide$Diet, sep="_")
R_median_wide <- dplyr::left_join(R_median_wide,R_pvalue, by=c("phylum_Diet"))

R_paintdata <- na.omit(R_median_wide)
R_paintdata$group <- ifelse(R_paintdata$pvalue<0.05, R_paintdata$aftertreat, "nosig")
R_paintdata$group <- factor(R_paintdata$group,levels = c("up","down","nosig")) #https://blog.csdn.net/weixin_56198196/article/details/124518774
##计算该饮食下百分比
paintdata <- R_paintdata %>%
  group_by(Diet) %>%
  count(group)%>%
  mutate(percent = n/sum(n))

##### Painting
ggbarplot(paintdata, x="group", y="percent",fill="group", 
               color=NA, #条边框无色
               palette = c("#F48892","#91CAE8", "grey"), 
               ggtheme= theme_test(),
               width = 0.6, #条宽
               position=position_dodge(.5), #条间距
               facet.by = "Diet", 
               ncol=8, nrow=1,
               short.panel.labs = T)
ggsave("../fig1-F-3.pdf") 
