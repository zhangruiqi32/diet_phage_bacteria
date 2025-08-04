##Script to Figure 1-B and Supplementary fig2-AB.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","ggradar","fmsb","ggplot2")
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
factorss2 <- read.csv("../data/fig1B.csv")

############################## Data processing ##############################
factorss2$level2 <- factor(factorss2$level,levels=c("VG","BG","BS","VS","BF","VF"),
                           labels=c("Viruses","Bacteria","Bacteria","Viruses","Bacteria","Viruses"))
factorss2$level3 <- factor(factorss2$level,levels=c("BG","BS","VG","VS","BF","VF"),
                           labels=c("Genus","Species","Genus","Species","Family","Family"))
factorss2 <- factorss2[order(factorss2$level2,factorss2$V4),]

############################## [Fig1-B] Painting #####################################
i=1
factorss3 <- factorss2[factorss2$Type==unique(factorss2$Type)[i],]
j=1 #j=1【BS VS】 j=2【BF VF】 j=3【VG BG】
factorss4 <- factorss3[factorss3$level3==unique(factorss3$level3)[j],]
factorss4 <- pivot_wider(factorss4[c("factor","V4","level")],names_from = "factor",values_from = "V4")
factorss4 <- column_to_rownames(factorss4, var = "level")

max_min <- data.frame(
  Sex= c(0.3, 0), Antibiotic = c(0.3, 0), Host = c(0.3, 0),
  Continent = c(0.3, 0), Country = c(0.3, 0), Diease = c(0.3, 0),
  Age2 = c(0.3, 0), Diet1 = c(0.3, 0), Nutrition1 = c(0.3, 0)
) #j=3【VG BG】时 max=0.35
rownames(max_min) <- c("Max", "Min")
df <- rbind(max_min, factorss4) 
colors <- c("#D7301F","#2171B5")
p<- radarchart(df,axistype = 1,
               # Customize the polygon
               pcol = colors,
               #pcol = colors,
               #pfcol = scales::alpha(colors, 0.2), #覆盖填充
               pty = 32, # 点的类型，16为默认实现圆点，32为无点
               plwd = 3,#线宽
               plty = 1,
               # Customize the grid
               cglcol = "grey20",
               cglty = 1,#线条类型 eg.lty=1表示实线，lty=2表示虚线，lty=3表示点线
               # Customize the axis
               axislabcol = "grey",
               caxislabels = seq(1,3,0.5),
               cglwd = 0.6,
               vlcex = 0.9,# 标签尺寸
)
legend(x = 1, y = 1.2, # legend的位置
       legend = rownames(df[3:4, ]), # legend为3—4行的行名称
       bty = "n",
       pch = 20 ,
       col = scales::alpha(colors, 0.6),
       text.col = "black",
       cex = 1,
       pt.cex = 2 )
ggsave("../fig1-B.pdf") 


############################## [FigS2-A] Painting #####################################
i=1
factorss3 <- factorss2[factorss2$Type==unique(factorss2$Type)[i],]
j=2 #j=1【BS VS】 j=2【BF VF】 j=3【VG BG】
factorss4 <- factorss3[factorss3$level3==unique(factorss3$level3)[j],]
factorss4 <- pivot_wider(factorss4[c("factor","V4","level")],names_from = "factor",values_from = "V4")
factorss4 <- column_to_rownames(factorss4, var = "level")

max_min <- data.frame(
  Sex= c(0.3, 0), Antibiotic = c(0.3, 0), Host = c(0.3, 0),
  Continent = c(0.3, 0), Country = c(0.3, 0), Diease = c(0.3, 0),
  Age2 = c(0.3, 0), Diet1 = c(0.3, 0), Nutrition1 = c(0.3, 0)
) #j=3【VG BG】时 max=0.35
rownames(max_min) <- c("Max", "Min")
df <- rbind(max_min, factorss4) 
colors <- c("#D7301F","#2171B5")
p<- radarchart(df,axistype = 1,
               # Customize the polygon
               pcol = colors,
               #pcol = colors,
               #pfcol = scales::alpha(colors, 0.2), #覆盖填充
               pty = 32, # 点的类型，16为默认实现圆点，32为无点
               plwd = 3,#线宽
               plty = 1,
               # Customize the grid
               cglcol = "grey20",
               cglty = 1,#线条类型 eg.lty=1表示实线，lty=2表示虚线，lty=3表示点线
               # Customize the axis
               axislabcol = "grey",
               caxislabels = seq(1,3,0.5),
               cglwd = 0.6,
               vlcex = 0.9,# 标签尺寸
)
legend(x = 1, y = 1.2, # legend的位置
       legend = rownames(df[3:4, ]), # legend为3—4行的行名称
       bty = "n",
       pch = 20 ,
       col = scales::alpha(colors, 0.6),
       text.col = "black",
       cex = 1,
       pt.cex = 2 )
ggsave("../figS2-A.pdf") 



############################## [FigS2-B] Painting #####################################
i=1
factorss3 <- factorss2[factorss2$Type==unique(factorss2$Type)[i],]
j=3 #j=1【BS VS】 j=2【BF VF】 j=3【VG BG】
factorss4 <- factorss3[factorss3$level3==unique(factorss3$level3)[j],]
factorss4 <- pivot_wider(factorss4[c("factor","V4","level")],names_from = "factor",values_from = "V4")
factorss4 <- column_to_rownames(factorss4, var = "level")

max_min <- data.frame(
  Sex= c(0.35, 0), Antibiotic = c(0.35, 0), Host = c(0.3, 0),
  Continent = c(0.35, 0), Country = c(0.35, 0), Diease = c(0.35, 0),
  Age2 = c(0.35, 0), Diet1 = c(0.35, 0), Nutrition1 = c(0.35, 0)
) #j=3【VG BG】时 max=0.35
rownames(max_min) <- c("Max", "Min")
df <- rbind(max_min, factorss4) 
colors <- c("#D7301F","#2171B5")
p<- radarchart(df,axistype = 1,
               # Customize the polygon
               pcol = colors,
               #pcol = colors,
               #pfcol = scales::alpha(colors, 0.2), #覆盖填充
               pty = 32, # 点的类型，16为默认实现圆点，32为无点
               plwd = 3,#线宽
               plty = 1,
               # Customize the grid
               cglcol = "grey20",
               cglty = 1,#线条类型 eg.lty=1表示实线，lty=2表示虚线，lty=3表示点线
               # Customize the axis
               axislabcol = "grey",
               caxislabels = seq(1,3,0.5),
               cglwd = 0.6,
               vlcex = 0.9,# 标签尺寸
)
legend(x = 1, y = 1.2, # legend的位置
       legend = rownames(df[3:4, ]), # legend为3—4行的行名称
       bty = "n",
       pch = 20 ,
       col = scales::alpha(colors, 0.6),
       text.col = "black",
       cex = 1,
       pt.cex = 2 )
ggsave("../figS2-B.pdf") 
