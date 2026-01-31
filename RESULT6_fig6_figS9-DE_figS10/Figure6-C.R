##Script to Figure 6-C.
############################## Packages loading ############################## 
library(ggplot2)

############################## Data loading ##############################
HGTs_bidui5<- read.table("~/cooperation/202409zhaofengxiang/data/fig6C_HGTs_bidui5")
reads2 <- read.table("~/cooperation/202409zhaofengxiang/data/fig6C_reads2.txt")

############################## Data processing ##############################
HGTs_bidui5 <- HGTs_bidui5[HGTs_bidui5$BioProject != "PRJNA704567" & HGTs_bidui5$Diet == "High_fat_diet", ]#去除PRJNA290729

reads2$Group <- as.character(reads2$Group)
reads2$Type <- ifelse(
  grepl("^HFD_D(10|14|18)$", reads2$Group), "control",
  ifelse(grepl("^FUC_D(10|14|18)$", reads2$Group), "case", NA)
)
reads2_filtered <- reads2[!is.na(reads2$Type) & reads2$Group %in% c("HFD_D10","HFD_D14","HFD_D18","FUC_D10","FUC_D14","FUC_D18"), ]

reads2_filtered$BioProject <- "Fucoidan"
reads2_filtered$Sample <- reads2_filtered$Sample_id
reads2_filtered <- reads2_filtered[reads2_filtered$Normalized_Mapped_reads!=0,]
reads2_filtered$Metabolism <- reads2_filtered$type
reads2_filtered$Nutrition1 <- reads2_filtered$Type
reads2_filtered$Mapped_number_nomalized2 <- log10(reads2_filtered$Normalized_Mapped_reads * 100)

HGTs_bidui5$Mapped_number_nomalized2 <- log10(HGTs_bidui5$Mapped_number_nomalized * 100)
common_cols <- c("BioProject", "Sample", "type", "Nutrition1", "Mapped_number_nomalized2")
HGT_combined <- rbind(HGTs_bidui5[, common_cols], reads2_filtered[, common_cols])
HGT_sampleAb <- aggregate(list(Mapped_number_nomalized = HGT_combined$Mapped_number_nomalized2),
                          by = list(BioProject = HGT_combined$BioProject,
                                    Sample = HGT_combined$Sample,
                                    Metabolism = HGT_combined$type,
                                    Nutrition1 = HGT_combined$Nutrition1),
                          sum)

HGT_sampleAb$Type <- as.character(factor(HGT_sampleAb$Nutrition1,
                                         levels = c("APS", "Baseline", "Control", "FOS", "Inulin", "case", "control"),
                                         labels = c("case", "control", "control", "case", "case", "case", "control")))

HGT_FC <- aggregate(list(mean_Ab = HGT_sampleAb$Mapped_number_nomalized),
                    by = list(BioProject = HGT_sampleAb$BioProject,
                              Metabolism = HGT_sampleAb$Metabolism,
                              Type = HGT_sampleAb$Type),
                    mean)

HGT_FCpaint <- tidyr::pivot_wider(HGT_FC, names_from = "Type", values_from = "mean_Ab")
HGT_FCpaint$log2FC <- log2(HGT_FCpaint$case / HGT_FCpaint$control)
HGT_FCpaint$FC <- HGT_FCpaint$case / HGT_FCpaint$control
HGT_FCpaint$BioProject <- as.factor(HGT_FCpaint$BioProject)
HGT_FCpaint <- dplyr::filter(HGT_FCpaint, Metabolism %in% c("integrase", "recombinase", "transposase"))

# 修改BioProject名称为全称
HGT_FCpaint$BioProject <- factor(HGT_FCpaint$BioProject,
                                 levels = c("PRJNA615253", "mgp6153", "PRJNA731974", "Fucoidan"),
                                 labels = c("Astragluspolysaccharide(APS)", 
                                            "Fructooligosaccharide(FOS)", 
                                            "Inulin", 
                                            "Fucoidan"))

############################## [Fig6-C] Painting #####################################
library(ggplot2)
ggplot(HGT_FCpaint, aes(x = Metabolism, y = log2FC, group = BioProject, color = BioProject)) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, size = 0.5, colour = "gray") +
  scale_color_manual(values = c("Astragluspolysaccharide(APS)" = "#2F5B60", 
                                "Fructooligosaccharide(FOS)" = "#9CB4AD", 
                                "Inulin" = "#AAD3C7", 
                                "Fucoidan" = "#5C887C")) +
  labs(y = "log2FC", x = "HGTs related genes", color = "Treatment") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())
