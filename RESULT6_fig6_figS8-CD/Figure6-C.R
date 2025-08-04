##Script to Figure 6-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggplot2")
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
HGTs_bidui5<- read.table("../data/fig6C_HGTs_bidui5")

############################## Data processing ##############################
HGTs_bidui5 <- HGTs_bidui5[HGTs_bidui5$Diet =="High_fat_diet",] #去除PRJNA290729

HGTs_bidui5$Mapped_number_nomalized2 <- log10(HGTs_bidui5$Mapped_number_nomalized*100)
HGT_sampleAb <- aggregate(list(Mapped_number_nomalized= HGTs_bidui5$Mapped_number_nomalized2 ),
                          by=list(BioProject=HGTs_bidui5$BioProject,
                                  Sample=HGTs_bidui5$Sample,
                                  Metabolism=HGTs_bidui5$type,
                                  Nutrition1=HGTs_bidui5$Nutrition1),
                          sum ) 
HGT_sampleAb$Type <- as.character(factor(HGT_sampleAb$Nutrition1,
                                         levels = c("APS", "Baseline" ,"Control",  "FOS", "Inulin","NR"),
                                         labels = c("case","control","control","case","case","case")))

HGT_FC <- aggregate(list(mean_Ab = HGT_sampleAb$Mapped_number_nomalized ),
                    by=list(BioProject=HGT_sampleAb$BioProject,
                            Metabolism=HGT_sampleAb$Metabolism,
                            Type=HGT_sampleAb$Type),
                    mean ) #等价于使用dplyr的group_by和summarise函数进行聚合求和
HGT_FCpaint <- tidyr::pivot_wider(HGT_FC,names_from = "Type",values_from = "mean_Ab")
HGT_FCpaint$log2FC <- log2( HGT_FCpaint$case / HGT_FCpaint$control )
HGT_FCpaint$FC <- HGT_FCpaint$case / HGT_FCpaint$control
HGT_FCpaint$BioProject <- as.factor(HGT_FCpaint$BioProject) 
HGT_FCpaint <- dplyr::filter(HGT_FCpaint, Metabolism %in% c("integrase",
                                                            "recombinase", 
                                                            "transposase"))


############################## [Fig6-C] Painting #####################################
ggplot(HGT_FCpaint,aes(x=Metabolism, y=log2FC, group = BioProject, color=BioProject))+
  geom_line()+
  geom_point(size=2)+
  geom_hline(yintercept=0, size=0.5, colour="gray")+ 
  scale_color_manual(values = c("PRJNA615253"="#2F5B60",
                                "PRJNA704567"="#5C887C",
                                "mgp6153"="#9CB4AD",
                                "PRJNA731974"="#AAD3C7"))+
  labs(y="log2FC", x="HGTs related genes")+
  theme_bw()+
  theme(panel.grid=element_blank())

ggsave("../fig6-C.pdf") 
