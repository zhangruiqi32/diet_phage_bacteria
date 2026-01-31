##Script to Supplementary fig6-D.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggplot2","ggalluvial")
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

library(ggalluvial)
library(ggplot2)
############################## Data loading ##############################
plot_data <- read.csv("~/cooperation/202409zhaofengxiang/data/figS6D_plot_data.csv",row.names = 1)

############################## Data processing ##############################
plot_data <- plot_data[plot_data$Count!=0,]
plot_data$Group <- factor(plot_data$Group,levels = c("CON","HFD","FUC"))
plot_data$Bact_p <- factor(plot_data$Bact_p,levels=rev(c("Bacteroidetes",
                                                         "Firmicutes",
                                                         "Actinobacteria","Proteobacteria","Cyanobacteria","Fusobacteria")))
############################## [FigS6-D] Painting #####################################
ggplot(plot_data,aes(x=Time,y=log10(Count), alluvium=Bact_p, stratum=Bact_p, fill=Bact_p))+
  geom_col(width=0.33)+
  geom_flow(aes(alluvium = Bact_p), 
            stat = "alluvium", 
            alpha = 0.2, 
            lode.guidance = "frontback")+
  facet_grid(~Group,scales="free_x",space="free")+
  scale_fill_manual(values=c("Bacteroidetes"="#85C3B9","Firmicutes"="#DDEAC4","Actinobacteria"="#DCDDDD","Proteobacteria"="#DCDDDD","Cyanobacteria"="#DCDDDD","Fusobacteria"="#DCDDDD"))+
  labs(x="Time point",y="log10(Number of Correlations)",fill="Phylum")+
  theme_test()
  

