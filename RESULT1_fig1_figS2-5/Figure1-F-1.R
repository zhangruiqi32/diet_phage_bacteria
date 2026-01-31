##Script to Figure 1-F-1.   

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("vegan","dplyr","grid","PupillometryR", "SuppDists", "ggpubr", "ggplot2")
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
shannon_paint <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F-1_shannon.csv", row.names = 1, header = T)


############################## [Fig1-F-1] Painting #####################################
ggviolin(shannon_paint, x="group", y="shannon",fill="group",
              palette = c("#2C7FB8","#D7301F"), alpha = 0.5,
              facet.by = "Diet", ncol=8, nrow=1,
              outlier.shape = NA,
              width = 0.5,
              add = "boxplot", 
              add.params = list(fill = "group", alpha=1, width = 0.1,linetype = 1)) + 
  stat_compare_means(label =  "p.signif", method="wilcox",
                            label.x.npc ="left", size = 5) 

ggplot(shannon_paint, aes(x=group, y=shannon))  +
  geom_flat_violin(aes(fill=group),position=position_nudge(x=.25),color=NA,alpha=0.5) +
  geom_jitter(size=1.2,aes(color=group),width=0.1,alpha=0.7,pch=16) + #pch=shape
  geom_boxplot(aes(fill=group),position=position_nudge(x=0.25),
               width=.1,size=0.5,outlier.color=NA) +
  #coord_flip() + #x和y翻转
  scale_fill_manual(values = c("#00b8e5","#FA8072"))+
  scale_color_manual(values = c("#00b8e5","#FA8072"))+
  theme_test() + #去除底纹网格线
  theme( axis.text = element_text(size=13),
         axis.title =  element_text(size=15),
         legend.position="none")+
  facet_wrap(~Diet,
                  ncol=8, nrow=1,
                  scales = "fixed")#fixed各子图坐标轴完全一样
ggsave("../fig1-F-1.pdf") 

