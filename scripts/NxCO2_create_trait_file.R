###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(LeafArea)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chlorophyll <- read.csv("../data_sheets/NxCO2_chlorophyllExtractions.csv")
co2.response.wk7 <- read.csv("../data_sheets/NxCO2_co2_resp_wk7.csv")
rd.wk7 <- read.csv("../data_sheets/NxCO2_rd_wk7.csv")

###############################################################################
## Load custom fxns
###############################################################################
source("/Users/eaperkowski/git/r_functions/temp_standardize.R")
source("/Users/eaperkowski/git/r_functions/calc_chi.R")
source("/Users/eaperkowski/git/r_functions/stomatal_limitation.R")

###############################################################################
## Fit curves
###############################################################################
aci.prep <- co2.response.wk7 %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair,	Tleaf) %>%
  arrange(id) %>%
  left_join(rd.wk7, by = "id") %>%
  dplyr::select(-week) %>%
  separate(col = "id",
           sep = "(_*)[_]_*",
           into = c("co2.trt", "inco", "n.trt", "rep"),
           remove = FALSE) %>%
  group_by(id) %>%
  mutate(rd25 = temp_standardize(rd,
                                 estimate.type = "Rd",
                                 pft = "C3H",
                                 standard.to = 25,
                                 tLeaf = Tleaf,
                                 tGrow = 22.5)) %>%
  data.frame()


## Calculate leaf disk area
ij.path <- "/Applications/ImageJ.app"
imagepath <- "/Users/eaperkowski/git/2022_NxCO2_growthChamber/leaf_area/chl_disk_scans/"

leaf.disk.area <- run.ij(path.imagej = ij.path,
                    set.directory = imagepath,
                    distance.pixel = 117.9034,
                    known.distance = 1,
                    set.memory = 300, low.size = 0.1)

names(leaf.disk.area)[1:2] <- c("id", "disk_area")

leaf.disk.area$id <- gsub("_disk", "", leaf.disk.area$id)

## Join leaf disk area with chlorophyll file
chlorophyll %>% full_join(leaf.disk.area) %>%
  mutate(chlA_ugml = (12.47 * a_665) - (3.62 * a_649),
         chlB_ugml = (25.06 * a_649) - (6.5 * a_665),
         chlA_mmol_m2 = chlA_ugml / 0.0001 / 8.935e8 * 1000,
         chlB_mmol_m2 = chlB_ugml / 0.0001 / 9.075e8 * 1000,
         chla.chlb = chlA_ugml / chlB_ugml) %>%
  arrange(chla.chlb)
## note mL converts to cm^3. Conversion factor: 1 cm^3 = 0.0001 m^2
## ChlA mol wgt = 893.5 g/mol = 8.935e8 μg/mol
## ChlB mol wgt = 907.5 g/mol = 9.075e8 μg/mol
## μg      1 mL       1 cm^3        1 mol ChlA      1000 mmol      mmol ChlA
## --- *  ------  * ----------- * --------------- * ---------- =  -----------
## mL     1 cm^3    0.0001 m^2    8.935e8 μg ChlA     1 mol           m^2
##
##
## μg      1 mL       1 cm^3        1 mol ChlB      1000 mmol      mmol ChlB
## --- *  ------  * ----------- * --------------- * ---------- =  -----------
## mL     1 cm^3    0.0001 m^2    9.075e8 μg ChlB     1 mol           m^2
