###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
co2.response.wk7 <- read.csv("../data_sheets/NxCO2_co2_resp_wk7.csv")
rd.wk7 <- read.csv("../data_sheets/NxCO2_rd_wk7.csv")
biomass_area <- read.csv("../data_sheets/NxCO2_tla_biomass_data.csv")




###############################################################################
## Prep data frame to fit A/Ci curves
###############################################################################
# Note: sub-setting to less than 700 μmol mol^-1 CO2 for LEMONTREE; will need 
# to 

aci.prep <- co2.response.wk7 %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair,	Tleaf) %>%
  arrange(id) %>%
  left_join(rd.wk7, by = "id") %>%
  dplyr::select(-week) %>%
  group_by(id) %>%
  mutate(rd25 = temp_standardize(rd,
                                 estimate.type = "Rd",
                                 pft = "C3H",
                                 standard.to = 25,
                                 tLeaf = Tleaf,
                                 tGrow = 22.5),
         keep.row = "yes") %>%
  data.frame()

aci.temps <- aci.prep %>%
  group_by(id) %>%
  summarize(tLeaf = mean(Tleaf, na.rm = TRUE))

###############################################################################
## Run A/Ci curves
###############################################################################

## Remove rows based on A/Ci fits, and also include all points measured
## at 0 ppm CO2. Workshop w/ Licor noted that 0ppm CO2 turns off mixing fan. 
#aci.merged$keep.row[aci.merged$A < -1.5] <- "no"
aci.prep$keep.row[c()] <- "no"

#####################################################################
# Coarse run through on A/Ci curves
#####################################################################
aci.fits <- aci.prep %>% filter(keep.row == "yes") %>%
  fitacis(group = "id",
          varnames = list(ALEAF = "A",
                          Tleaf = "Tleaf",
                          Ci = "Ci",
                          PPFD = "Qin", 
                          Rd = "rd25"),
          fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

summary(aci.fits)

photo.params <- coef(aci.fits) %>%
  select(id:Rd) %>%
  full_join(aci.temps) %>%
  mutate(vcmax25 = temp_standardize(Vcmax, "Vcmax", standard.to = 25,
                                    tLeaf = tLeaf, tGrow = 22.5),
         jmax25 = temp_standardize(Jmax, "Jmax", standard.to = 25,
                                   tLeaf = tLeaf, tGrow = 22.5),
         jmax25 = ifelse(jmax25 < 0, NA, jmax25),
         jmax25.vcmax25 = jmax25 / vcmax25,) %>%
  slice(-c(25, 51, 59, 91))

## Change incorrect IDs
photo.params$id[photo.params$id == "e_n_280_26"] <- "e_y_280_26"
photo.params$id[photo.params$id == "e_y_350_68"] <- "e_n_350_68" 
photo.params$id[photo.params$id == "a_y_280_133"] <- "a_n_280_133" 
photo.params$id[photo.params$id == "a_y_105_123"] <- "a_n_105_123" 
photo.params$id[photo.params$id == "a_y_70_79"] <- "a_y_35_79" 