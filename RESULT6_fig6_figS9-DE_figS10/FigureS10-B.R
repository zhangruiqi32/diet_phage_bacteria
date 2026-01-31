##Script to Supplementary fig10-B.

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
plot_data <- read.csv("~/cooperation/202409zhaofengxiang/data/figS10B_plot_data.csv",row.names = 1)

############################## Data processing ##############################
#plot_data <- plot_data[plot_data$Count!=0,]
plot_data$Bioproject <- factor(plot_data$Bioproject,levels = c("PRJNA615253","PRJNA731974","mgp6153"))
plot_data$Treat <- factor(plot_data$Treat,levels = c("Control","APS","Inulin","FOS"))
plot_data$Bact_p <- factor(plot_data$Bact_p,levels=rev(c("Bacteroidetes",
                                                         "Firmicutes",
                                                         "Actinobacteria","Proteobacteria","Cyanobacteria","Deinococcus-Thermus")))
plot_data <- arrange(plot_data,Bioproject,Treat,Bact_p)
############################## [FigS10-B] Painting #####################################
ggplot(plot_data,aes(x=Treat,y=Count, fill=Bact_p)) + 
         geom_col(width=0.33)+ 
  geom_flow(aes(alluvium = Bact_p), 
            stat = "alluvium", 
            alpha = 0.2, 
            lode.guidance = "frontback")+
         facet_wrap(~Bioproject, scales="free") + 
  scale_fill_manual(values=c("Bacteroidetes"="#85C3B9","Firmicutes"="#DDEAC4","Actinobacteria"="#DCDDDD","Proteobacteria"="#DCDDDD","Cyanobacteria"="#DCDDDD","Deinococcus-Thermus"="#DCDDDD"))+
         labs(x="",y="Number of Correlations",fill="Phylum")+
         theme_test()

