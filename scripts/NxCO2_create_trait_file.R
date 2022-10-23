###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(LeafArea)
library(plantecophys)
library(car)
library(emmeans)

###############################################################################
## Import files
###############################################################################
biomass_area <- read.csv("../data_sheets/NxCO2_tla_biomass_data.csv")
ids <- data.frame(id = biomass_area$id)
ids <- separate(ids, id, sep = "_", into = c("co2", "inoc", "n.trt", "rep"),
                remove = FALSE)

###############################################################################
## Chlorophyll content
###############################################################################
## Load chlorophyll plate reader file
chlor.files <- list.files(path = "../data_sheets/chlorophyll/", 
                          pattern = "\\.txt$", recursive = TRUE,
                          full.names = TRUE)
chlor.files <- setNames(chlor.files, chlor.files)

chlor.df <- lapply(chlor.files, read.delim) %>% 
  reshape::merge_all() %>%
  filter(id != "blank") %>%
  select(rep = id, everything()) %>%
  group_by(rep) %>%
  summarize(abs.649 = mean(abs_649_corrected, na.rm = TRUE),
            cv.649 = (sd(abs_649_corrected, na.rm = TRUE) / abs.649) * 100,
            abs.665 = mean(abs_665_corrected, na.rm = TRUE),
            cv.665 = (sd(abs_665_corrected, na.rm = TRUE) / abs.665) * 100) %>%
  mutate(abs.649 = ifelse(abs.649 < 0, 0, abs.649),
         abs.665 = ifelse(abs.665 < 0, 0, abs.665))
chlor.df <- filter(chlor.df, rep != "12")
chlor.df$rep[chlor.df$rep == "12_real"] <- "12"

# Replace rep number with full plant ID
chlor.df <- chlor.df %>%
  full_join(ids) %>%
  select(id, abs.649:cv.665)

## Calculate leaf disk area
ij.path <- "/Applications/ImageJ.app"
imagepath.disk <- "/Users/eaperkowski/git/2022_NxCO2xI/leaf_area/chl_disk_scans/"
imagepath.chlLeaf <- "/Users/eaperkowski/git/2022_NxCO2xI/leaf_area/chl_leaf_scans/"

chlor.disk.area <- run.ij(path.imagej = ij.path,
                         set.directory = imagepath.disk,
                         distance.pixel = 117.9034,
                         known.distance = 1, low.size = 0.05)
chlor.leaf.area <- run.ij(path.imagej = ij.path,
                          set.directory = imagepath.chlLeaf,
                          distance.pixel = 117.9034,
                          known.distance = 1, low.size = 0.05)

## Rename area column, merge into central file
names(chlor.disk.area)[1:2] <- c("id", "disk.area")
names(chlor.leaf.area)[1:2] <- c("id", "chl.leaf.area")

chlor.disk.area$id <- gsub("_disk", "", chlor.disk.area$id)
chlor.leaf.area$id <- gsub("_chlLeaf", "", chlor.leaf.area$id)

chl.leaf.area <- chlor.disk.area %>%
  full_join(chlor.leaf.area)

## Join leaf disk area with chlorophyll file
chlorophyll <- chlor.df %>% 
  full_join(chl.leaf.area, by = "id") %>%
  full_join(biomass_area) %>%
  select(id, abs.649, abs.665, chl.biomass = chlor_biomass, disk.area, chl.leaf.area) %>%
  mutate(chlA_ugml = (12.19 * abs.665) - (3.56 * abs.649),
         chlB_ugml = (21.99 * abs.649) - (5.32 * abs.665),
         chlA_ugml = ifelse(chlA_ugml < 0, 0, chlA_ugml),
         chlB_ugml = ifelse(chlB_ugml < 0, 0 , chlB_ugml),
         chlA_gml = chlA_ugml / 1000000,
         chlB_gml = chlB_ugml / 1000000,
         chlA_g = chlA_gml * 10, # extracted in 10mL DMSO
         chlB_g = chlB_gml * 10, # extracted in 10mL DMSO
         chlA_gm2 = chlA_g / (disk.area / 10000),
         chlB_gm2 = chlB_g / (disk.area / 10000),
         chlA_mmolm2 = chlA_gm2 / 893.51 * 1000,
         chlB_mmolm2 = chlB_gm2 / 907.47 * 1000,
         chl_marea = chl.biomass / (chl.leaf.area / 10000),
         chlA_mmolg = chlA_mmolm2 / chl_marea,
         chlB_mmolg = chlB_mmolm2 / chl_marea,
         chl_mmolg = chlA_mmolg + chlB_mmolg,
         chl_mmolm2 = chlA_mmolm2 + chlB_mmolm2,
         chlA.chlB = chlA_g / chlB_g) %>%
  select(id, chlA_ugml:chlA.chlB)

###############################################################################
## Focal leaf area, leaf nitrogen content
###############################################################################
## Load CN data files
cn.files <- list.files(path = "../data_sheets/costech_results/", 
                          pattern = "\\.csv$", recursive = TRUE,
                          full.names = TRUE)
cn.files <- setNames(cn.files, cn.files)

focal.nmass <- lapply(cn.files, read.csv) %>% reshape::merge_all() %>%
  filter(type == "unknown") %>%
  mutate(id = gsub("_focal", "", id)) %>%
  select(id, nmass = n.weight.percent)

## Calculate leaf disk area
imagepath.focal <- "/Users/eaperkowski/git/2022_NxCO2xI/leaf_area/focal_scans/"

focal.area <- run.ij(path.imagej = ij.path,
                     set.directory = imagepath.focal,
                     distance.pixel = 117.9034,
                     known.distance = 1, low.size = 0.05,
                     set.memory = 30)

names(focal.area)[1:2] <- c("id", "focal.area")
focal.area$id <- gsub("_focal", "", focal.area$id)

leafn <- focal.area %>% full_join(biomass_area) %>%
  full_join(focal.nmass) %>%
  select(id, focal.area, focal.biomass = focal_biomass, nmass) %>%
  mutate(marea = focal.biomass / (focal.area / 10000),
         narea = (nmass/100) * marea)
