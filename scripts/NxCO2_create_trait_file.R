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
chlorophyll <- read.csv("../data_sheets/NxCO2_chlorophyllExtractions.csv")
co2.response.wk7 <- read.csv("../data_sheets/NxCO2_co2_resp_wk7.csv")
rd.wk7 <- read.csv("../data_sheets/NxCO2_rd_wk7.csv")
biomass_area <- read.csv("../data_sheets/NxCO2_tla_biomass_data.csv")

###############################################################################
## Load custom fxns
###############################################################################
source("/Users/eaperkowski/git/r_functions/temp_standardize.R")
source("/Users/eaperkowski/git/r_functions/calc_chi.R")
source("/Users/eaperkowski/git/r_functions/stomatal_limitation.R")

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
  filter(Ci < 600) %>%
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
          fitTPU = FALSE, Tcorrect = FALSE, useRd = FALSE)

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

photo.params <- photo.params %>%
  separate(id, into = c("co2", "inoc", "n.trt", "rep")) %>%
  unite(col = "id", co2:rep, remove = FALSE) %>%
  full_join(biomass_area, by = "id")

write.csv(photo.params, "../data_sheets/NxCO2_datasheet.csv", row.names = FALSE)

## Preliminary models for LEMONTREE meeting
photo.params$vcmax25[91] <- NA

vcmax25 <- lm(vcmax25 ~ as.factor(co2) * inoc * as.numeric(n.trt),
              data = photo.params)
summary(vcmax25)
Anova(vcmax25)

shapiro.test(residuals(vcmax25))
outlierTest(vcmax25)

## Preliminary plots for LEMONTREE meeting
inoc.labs <- c("Not inoculated", "Inoculated")
names(inoc.labs) <- c("n", "y")


ggplot(data = photo.params, aes(x = as.numeric(n.trt), y = vcmax25, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 30)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("V")["cmax25"]*" (μmol m"^"-2"*"s"^"-1"*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  theme_bw(base_size = 18)

ggplot(data = photo.params, aes(x = as.numeric(n.trt), y = jmax25, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = "lm") +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("J")["max25"]*" (μmol m"^"-2"*"s"^"-1"*")")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  theme_bw(base_size = 18)


ggplot(data = photo.params, aes(x = as.numeric(n.trt), y = jmax25.vcmax25, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = "lm") +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  scale_y_continuous(limits = c(1.4,2), breaks = seq(1.4,2,0.2)) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("J")["max25"]*" (μmol m"^"-2"*"s"^"-1"*")")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  theme_bw(base_size = 18)



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
