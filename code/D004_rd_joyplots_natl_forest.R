# Load the packages
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(ggridges)

# Set working directory
path = c("/Users/zpzuo/My Drive/asrl_output/", "G://My Drive//asrl_output//", "/projectnb/amazondr/data10/cliveg/zpzuo/paper2_intermediate_data/")
wd = path[3]
setwd(wd)

# Load tree height raster
## Function
rasterToDf = function(file_path){
  img = raster(file_path)
  img_df = as.data.frame(img) %>%
    na.omit() 
  return(img_df)
}

# Relative Deviation (RD) - pCTH vs GEDI L3
## RD with GEDI -- FTG 100, national forest
rd_wGEDI_ftg100 = rasterToDf(paste0(wd, "RelDev_FTG100_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'White/Red/Jack Pine', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg100) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 120, national forest
rd_wGEDI_ftg120 = rasterToDf(paste0(wd, "RelDev_FTG120_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Spruce/Fir', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg120) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 140, national forest
rd_wGEDI_ftg140 = rasterToDf(paste0(wd, "RelDev_FTG140_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Longleaf/Slash Pine', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg140) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 160, national forest
rd_wGEDI_ftg160 = rasterToDf(paste0(wd, "RelDev_FTG160_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Loblolly/Shortleaf Pine', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg160) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 400, national forest
rd_wGEDI_ftg400 = rasterToDf(paste0(wd, "RelDev_FTG400_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Pine', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg400) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 500, national forest
rd_wGEDI_ftg500 = rasterToDf(paste0(wd, "RelDev_FTG500_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Hickory', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg500) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 600, national forest
rd_wGEDI_ftg600 = rasterToDf(paste0(wd, "RelDev_FTG600_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Gum/Cypress', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg600) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 800, national forest
rd_wGEDI_ftg800 = rasterToDf(paste0(wd, "RelDev_FTG800_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Maple/Beech/Birch', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg800) = c('RD', 'FTG', 'comparison')
## RD with GEDI -- FTG 900, national forest
rd_wGEDI_ftg900 = rasterToDf(paste0(wd, "RelDev_FTG900_wGEDI_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Aspen/Birch', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wGEDI_ftg900) = c('RD', 'FTG', 'comparison')


# Relative Deviation (RD) - pCTH vs Lang
## RD with Lang -- FTG 100, national forest
rd_wLang_ftg100 = rasterToDf(paste0(wd, "RelDev_FTG100_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'White/Red/Jack Pine', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wLang_ftg100) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 120, national forest
rd_wLang_ftg120 = rasterToDf(paste0(wd, "RelDev_FTG120_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Spruce/Fir', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg120) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 140, national forest
rd_wLang_ftg140 = rasterToDf(paste0(wd, "RelDev_FTG140_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Longleaf/Slash Pine', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg140) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 160, national forest
rd_wLang_ftg160 = rasterToDf(paste0(wd, "RelDev_FTG160_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Loblolly/Shortleaf Pine', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg160) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 400, national forest
rd_wLang_ftg400 = rasterToDf(paste0(wd, "RelDev_FTG400_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Pine', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg400) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 500, national forest
rd_wLang_ftg500 = rasterToDf(paste0(wd, "RelDev_FTG500_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Hickory', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg500) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 600, national forest
rd_wLang_ftg600 = rasterToDf(paste0(wd, "RelDev_FTG600_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Gum/Cypress', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg600) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 800, national forest
rd_wLang_ftg800 = rasterToDf(paste0(wd, "RelDev_FTG800_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Maple/Beech/Birch', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg800) = c('RD', 'FTG', 'comparison')
## Rd with Lang -- FTG 900, national forest
rd_wLang_ftg900 = rasterToDf(paste0(wd, "RelDev_FTG900_wLang_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Aspen/Birch', comparison = '(pCTH-Lang)/Lang')
colnames(rd_wLang_ftg900) = c('RD', 'FTG', 'comparison')


# Relative Deviation (RD) - pCTH vs Tolan
## RD with Tolan -- FTG 100, national forest
rd_wTolan_ftg100 = rasterToDf(paste0(wd, "RelDev_FTG100_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'White/Red/Jack Pine', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg100) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 120, national forest
rd_wTolan_ftg120 = rasterToDf(paste0(wd, "RelDev_FTG120_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Spruce/Fir', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg120) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 140, national forest
rd_wTolan_ftg140 = rasterToDf(paste0(wd, "RelDev_FTG140_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Longleaf/Slash Pine', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg140) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 160, national forest
rd_wTolan_ftg160 = rasterToDf(paste0(wd, "RelDev_FTG160_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Loblolly/Shortleaf Pine', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg160) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 400, national forest
rd_wTolan_ftg400 = rasterToDf(paste0(wd, "RelDev_FTG400_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Pine', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg400) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 500, national forest
rd_wTolan_ftg500 = rasterToDf(paste0(wd, "RelDev_FTG500_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Hickory', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg500) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 600, national forest
rd_wTolan_ftg600 = rasterToDf(paste0(wd, "RelDev_FTG600_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Gum/Cypress', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg600) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 800, national forest
rd_wTolan_ftg800 = rasterToDf(paste0(wd, "RelDev_FTG800_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Maple/Beech/Birch', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg800) = c('RD', 'FTG', 'comparison')
## RD with Tolan -- FTG 900, national forest
rd_wTolan_ftg900 = rasterToDf(paste0(wd, "RelDev_FTG900_wTolan_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Aspen/Birch', comparison = '(pCTH-Tolan)/Tolan')
colnames(rd_wTolan_ftg900) = c('RD', 'FTG', 'comparison')


# Relative Deviation (RD) - pCTH vs TreeMap
## RD with TreeMap -- FTG 100, national forest
rd_wTreeMap_ftg100 = rasterToDf(paste0(wd, "RelDev_FTG100_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'White/Red/Jack Pine', comparison = '(pCTH-GEDI)/GEDI')
colnames(rd_wTreeMap_ftg100) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 120, national forest
rd_wTreeMap_ftg120 = rasterToDf(paste0(wd, "RelDev_FTG120_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Spruce/Fir', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg120) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 140, national forest
rd_wTreeMap_ftg140 = rasterToDf(paste0(wd, "RelDev_FTG140_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Longleaf/Slash Pine', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg140) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 160, national forest
rd_wTreeMap_ftg160 = rasterToDf(paste0(wd, "RelDev_FTG160_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Loblolly/Shortleaf Pine', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg160) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 400, national forest
rd_wTreeMap_ftg400 = rasterToDf(paste0(wd, "RelDev_FTG400_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Pine', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg400) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 500, national forest
rd_wTreeMap_ftg500 = rasterToDf(paste0(wd, "RelDev_FTG500_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Hickory', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg500) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 600, national forest
rd_wTreeMap_ftg600 = rasterToDf(paste0(wd, "RelDev_FTG600_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Gum/Cypress', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg600) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 800, national forest
rd_wTreeMap_ftg800 = rasterToDf(paste0(wd, "RelDev_FTG800_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Maple/Beech/Birch', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg800) = c('RD', 'FTG', 'comparison')
## Rd with TreeMap -- FTG 900, national forest
rd_wTreeMap_ftg900 = rasterToDf(paste0(wd, "RelDev_FTG900_wTreeMap_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Aspen/Birch', comparison = '(pCTH-TreeMap)/TreeMap')
colnames(rd_wTreeMap_ftg900) = c('RD', 'FTG', 'comparison')

# Relative Deviation (RD) - pCTH vs Potapov
## RD with Potapov -- FTG 100, national forest
rd_wPotapov_ftg100 = rasterToDf(paste0(wd, "RelDev_FTG100_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'White/Red/Jack Pine', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg100) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 120, national forest
rd_wPotapov_ftg120 = rasterToDf(paste0(wd, "RelDev_FTG120_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Spruce/Fir', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg120) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 140, national forest
rd_wPotapov_ftg140 = rasterToDf(paste0(wd, "RelDev_FTG140_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Longleaf/Slash Pine', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg140) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 160, national forest
rd_wPotapov_ftg160 = rasterToDf(paste0(wd, "RelDev_FTG160_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Loblolly/Shortleaf Pine', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg160) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 400, national forest
rd_wPotapov_ftg400 = rasterToDf(paste0(wd, "RelDev_FTG400_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Pine', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg400) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 500, national forest
rd_wPotapov_ftg500 = rasterToDf(paste0(wd, "RelDev_FTG500_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Hickory', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg500) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 600, national forest
rd_wPotapov_ftg600 = rasterToDf(paste0(wd, "RelDev_FTG600_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Oak/Gum/Cypress', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg600) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 800, national forest
rd_wPotapov_ftg800 = rasterToDf(paste0(wd, "RelDev_FTG800_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Maple/Beech/Birch', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg800) = c('RD', 'FTG', 'comparison')
## Rd with Potapov -- FTG 900, national forest
rd_wPotapov_ftg900 = rasterToDf(paste0(wd, "RelDev_FTG900_wPotapov_maskS1NatlFor_1km.tif")) %>%
  mutate(FTG = 'Aspen/Birch', comparison = '(pCTH-Potapov)/Potapov')
colnames(rd_wPotapov_ftg900) = c('RD', 'FTG', 'comparison')

# Assemble RD rasters to a dataframe: 
# Column 1 - RD values
# Column 2 - Forest type group (FTG)
# Column 2 - comparison pair (PCTH vs GEDI, or PCTH vs Lang)
rd_wGEDI = bind_rows(rd_wGEDI_ftg100, rd_wGEDI_ftg120, rd_wGEDI_ftg140, rd_wGEDI_ftg160, rd_wGEDI_ftg400, rd_wGEDI_ftg500, rd_wGEDI_ftg600, rd_wGEDI_ftg800, rd_wGEDI_ftg900)
rd_wLang = bind_rows(rd_wLang_ftg100, rd_wLang_ftg120, rd_wLang_ftg140, rd_wLang_ftg160, rd_wLang_ftg400, rd_wLang_ftg500, rd_wLang_ftg600, rd_wLang_ftg800, rd_wLang_ftg900)
rd_wTolan = bind_rows(rd_wTolan_ftg100, rd_wTolan_ftg120, rd_wTolan_ftg140, rd_wTolan_ftg160, rd_wTolan_ftg400, rd_wTolan_ftg500, rd_wTolan_ftg600, rd_wTolan_ftg800, rd_wTolan_ftg900)
rd_wTreeMap = bind_rows(rd_wTreeMap_ftg100, rd_wTreeMap_ftg120, rd_wTreeMap_ftg140, rd_wTreeMap_ftg160, rd_wTreeMap_ftg400, rd_wTreeMap_ftg500, rd_wTreeMap_ftg600, rd_wTreeMap_ftg800, rd_wTreeMap_ftg900)
rd_wPotapov = bind_rows(rd_wPotapov_ftg100, rd_wPotapov_ftg120, rd_wPotapov_ftg140, rd_wPotapov_ftg160, rd_wPotapov_ftg400, rd_wPotapov_ftg500, rd_wPotapov_ftg600, rd_wPotapov_ftg800, rd_wPotapov_ftg900)


# Draw ridgeline plots
ridge_rd_wGEDI = ggplot(rd_wGEDI, aes(x=RD, y=FTG, fill=FTG)) +
  geom_density_ridges(alpha=0.6) +
  scale_fill_manual(values=c("White/Red/Jack Pine"="#a80084", "Spruce/Fir"="#ff73df", "Longleaf/Slash Pine"="#a80000", "Loblolly/Shortleaf Pine"="#ffaa00", "Oak/Pine"="#70a800", 
                             "Oak/Hickory"="#aaff00", "Oak/Gum/Cypress"="#73ffdf", "Maple/Beech/Birch"="#004da8", "Aspen/Birch"="#828282")) +
  theme_classic(base_size = 22) +
  labs(x='Relative Deviation: (pCTH-GEDI)/GEDI', y='') +
  coord_cartesian(xlim=c(-2.5,2.5)) + 
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) 
ridge_rd_wGEDI


ridge_rd_wLang = ggplot(rd_wLang[rd_wLang[,1]<5,], aes(x=RD, y=FTG, fill=FTG)) +
  geom_density_ridges(alpha=0.6) +
  scale_fill_manual(values=c("White/Red/Jack Pine"="#a80084", "Spruce/Fir"="#ff73df", "Longleaf/Slash Pine"="#a80000", "Loblolly/Shortleaf Pine"="#ffaa00", "Oak/Pine"="#70a800", 
                             "Oak/Hickory"="#aaff00", "Oak/Gum/Cypress"="#73ffdf", "Maple/Beech/Birch"="#004da8", "Aspen/Birch"="#828282")) +
  theme_classic(base_size = 22) +
  labs(x='Relative Deviation: (pCTH-Lang)/Lang', y='') +
  coord_cartesian(xlim=c(-2.5,2.5)) + 
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) 
ridge_rd_wLang


ridge_rd_wTolan = ggplot(rd_wTolan, aes(x=RD, y=FTG, fill=FTG)) +
  geom_density_ridges(alpha=0.6) +
  scale_fill_manual(values=c("White/Red/Jack Pine"="#a80084", "Spruce/Fir"="#ff73df", "Longleaf/Slash Pine"="#a80000", "Loblolly/Shortleaf Pine"="#ffaa00", "Oak/Pine"="#70a800", 
                             "Oak/Hickory"="#aaff00", "Oak/Gum/Cypress"="#73ffdf", "Maple/Beech/Birch"="#004da8", "Aspen/Birch"="#828282")) +
  theme_classic(base_size = 22) +
  labs(x='Relative Deviation: (pCTH-Tolan)/Tolan', y='') +
  coord_cartesian(xlim=c(-2.5,2.5)) + 
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) 
ridge_rd_wTolan


ridge_rd_wTreeMap = ggplot(rd_wTreeMap, aes(x=RD, y=FTG, fill=FTG)) +
  geom_density_ridges(alpha=0.6) +
  scale_fill_manual(values=c("White/Red/Jack Pine"="#a80084", "Spruce/Fir"="#ff73df", "Longleaf/Slash Pine"="#a80000", "Loblolly/Shortleaf Pine"="#ffaa00", "Oak/Pine"="#70a800", 
                             "Oak/Hickory"="#aaff00", "Oak/Gum/Cypress"="#73ffdf", "Maple/Beech/Birch"="#004da8", "Aspen/Birch"="#828282")) +
  theme_classic(base_size = 22) +
  labs(x='Relative Deviation: (pCTH-TreeMap)/TreeMap', y='') +
  coord_cartesian(xlim=c(-2.5,2.5)) + 
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) 
ridge_rd_wTreeMap


ridge_rd_wPotapov = ggplot(rd_wPotapov, aes(x=RD, y=FTG, fill=FTG)) +
  geom_density_ridges(alpha=0.6) +
  scale_fill_manual(values=c("White/Red/Jack Pine"="#a80084", "Spruce/Fir"="#ff73df", "Longleaf/Slash Pine"="#a80000", "Loblolly/Shortleaf Pine"="#ffaa00", "Oak/Pine"="#70a800", 
                             "Oak/Hickory"="#aaff00", "Oak/Gum/Cypress"="#73ffdf", "Maple/Beech/Birch"="#004da8", "Aspen/Birch"="#828282")) +
  theme_classic(base_size = 22) +
  labs(x='Relative Deviation: (pCTH-Potapov)/Potapov', y='') +
  coord_cartesian(xlim=c(-2.5,2.5)) + 
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) 
ridge_rd_wPotapov
