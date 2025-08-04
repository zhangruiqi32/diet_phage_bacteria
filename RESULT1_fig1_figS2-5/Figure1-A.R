##Script to Figure 1-A.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","sf","rnaturalearth", "rnaturalearthdata", "ggplot2")
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
df <- read_csv("../data/fig1A_map.csv") #sample
df_diet <- read_csv("../data/fig1A_mapdiet.csv") #diet label
df_nutrition <- read_csv("../data/fig1A_mapnutrition.csv") # component label

############################## Data processing ##############################
## 地图数据
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  select(subunit,continent,geometry)
## 样本数据
map <- df %>% 
  group_by(country) %>% 
  summarise(count=n(), sum(N)) 
names(map) <- c("country","project_N","sample_N")
## 合并成大洲数据
map2 <- left_join(
  map,world,
  by=c("country"="subunit")
)
map2 <- map2[,-5]
map_continent <- map2 %>% 
  group_by(continent) %>% 
  summarise(continent_project=sum(project_N),continent_sample=sum(sample_N)) 
#合并绘图数据
map3 <- left_join(
  world,
  map_continent,
  by=c("continent"="continent")
)
#测试颜色
map3 %>%
  ggplot() +
  geom_sf(aes(fill=continent_project, geometry=geometry), color=NA, size=0) +
  scale_fill_gradient(name = "# of Bioproject",
                      low = "#DEEBF7",
                      high = "#4292C6",
                      na.value="#F0F0F0",
                      limits=c(0,31),
                      breaks=seq(0,31,6))+
  theme_void() +
  theme(panel.grid = element_blank()) 
##在每块国家上加个文字
world_points<- st_centroid(world) #使用包 sf 中的 st_centroid 定义每块的质心
world_XY <- data.frame(world$subunit,st_coordinates(st_centroid(world$geometry)))
map_name <- left_join(
  map, #这里不能用替换UK和USA的 需要用原始map
  world_XY,
  by=c("country"="world.subunit")
) 
map_name$lable <- paste("(n = ",map_name$sample_N,")",sep = "")
#在每个名字上面标注饮食+营养种类数量
 
map_diet <- left_join(
  df_diet,
  world_XY,
  by=c("country"="world.subunit")
) 

map_nutrition <- left_join(
  df_nutrition,
  world_XY,
  by=c("country"="world.subunit")
) 

write.csv(map_diet,"../map_diet.csv") #change position of diet label
map_diet_XY <- read_csv("../data/fig1A_mapdiet_XY.csv") 
write.csv(map_nutrition,"../map_nutrition.csv") #change position of component label
map_nutrition_XY <- read_csv("../data/fig1A_mapnutrition_XY.csv") 


############################## [Fig1-A] Painting #####################################
#名字在圆点上面 样本数在圆点下面 饮食种类在n附近
ggplot() +
  geom_sf(data= map3,
          aes(fill=continent_project, geometry=geometry), color=NA, size=0) +
  scale_fill_gradient(name = "# of Bioproject",
                      low = "#BDD1E4",
                      high = "#548CB6",
                      na.value="#F0F0F0",
                      limits=c(0,26),
                      breaks=seq(0,26,6))+
  geom_point(data= map_name,aes(x=X, y=Y, label=country ),
             shape = 19, size = 0.8, color = "#8B0000") +
  geom_text(data= map_name,aes(x=X, y=Y+2.2, label=country),
            size = 2,
            color = "grey20", fontface = "plain", check_overlap = FALSE) +
  geom_text(data= map_name,aes(x=X, y=Y-1.8, label=lable),
            size = 1.5,
            color = "grey20", fontface = "plain", check_overlap = FALSE) +
  geom_point(data= map_diet_XY,aes(x=X, y=Y+4.1, color = score ),
             shape = 18, size = 0.9) +  
  scale_color_gradient(name = "# of DietType",
                       low =  "#FDBB84",
                       high =   "#D7301F",
                       limits=c(0,15),
                       breaks=seq(0,15,5))+
  theme_void() +
  theme(panel.grid = element_blank()) 
ggsave("../fig1A_map_diet.pdf", width = 12, height = 6) 
#Scale for colour is already present.
#Adding another scale for colour, which will replace the existing scale.
ggplot() +
  geom_sf(data= map3,
          aes(fill=continent_project, geometry=geometry), color=NA, size=0) +
  scale_fill_gradient(name = "# of Bioproject",
                      low = "#BDD1E4",
                      high = "#548CB6",
                      na.value="#F0F0F0",
                      limits=c(0,26),
                      breaks=seq(0,26,6))+
  geom_point(data= map_name,aes(x=X, y=Y, label=country ),
             shape = 19, size = 0.8, color = "#8B0000") +
  geom_text(data= map_name,aes(x=X, y=Y+2.2, label=country),
            size = 2,
            color = "grey20", fontface = "plain", check_overlap = FALSE) +
  geom_text(data= map_name,aes(x=X, y=Y-1.8, label=lable),
            size = 1.5,
            color = "grey20", fontface = "plain", check_overlap = FALSE) +
  geom_point(data= map_nutrition_XY,aes(x=X, y=Y-3.2, color = score ),
             shape = 18, size = 0.9) +  
  scale_color_gradient(name = "# of DietComponent",
                       low =  "#ADDD8E",
                       high =   "#006837",
                       limits=c(0,20),
                       breaks=seq(0,20,5))+
  theme_void() +
  theme(panel.grid = element_blank()) 
ggsave("../fig1A_map_nutrition.pdf", width = 12, height = 6) 

