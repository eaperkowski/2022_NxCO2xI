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
id <- read.csv("../data_sheets/NxCO2_id_datasheet.csv")
photo <- read.csv("../data_sheets/NxCO2_photo_data.csv")
photo.wk6 <- read.csv("../data_sheets/NxCO2_photo_data_wk6.csv")
isotopes <- read.csv("../data_sheets/NxCO2xI_isotope_data.csv")
d13c.air <- read.csv("../data_sheets/NxCO2xI_d13c_air.csv")


###############################################################################
## Load propN functions
###############################################################################
source("../../r_functions/propN_funcs.R")
source("../../r_functions/stomatal_limitation.R")
source("../../r_functions/calc_chi.R")
source("../../r_functions/calc_beta.R")
source("../../r_functions/calc_ndfa.R")

###############################################################################
## Calculate d13C in air between the two CO2 treatments
###############################################################################
d13c_air <- d13c.air %>%
  separate(id, into = c("name", "chamber", "co2", "rep")) %>%
  group_by(co2) %>%
  summarize(d13c.air = mean(d13c_air),
            co2.ppm = mean(co2_ppm)) %>%
  select(co2_cat = co2, everything()) %>%
  data.frame()

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
  dplyr::select(rep = id, everything()) %>%
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
chlor.df <- id %>%
  separate(id, c("co2", "inoc", "n.trt", "rep"), remove = FALSE) %>%
  full_join(chlor.df) %>%
  dplyr::select(id, abs.649, abs.665)

## Calculate leaf disk area
ij.path <- "/Applications/ImageJ.app"
imagepath.disk <- "/Users/eaperkowski/git/2022_NxCO2xI/leaf_area/chl_disk_scans/"
imagepath.chlLeaf <- "/Users/eaperkowski/git/2022_NxCO2xI/leaf_area/chl_leaf_scans/"

chlor.disk.area <- run.ij(path.imagej = ij.path,
                         set.directory = imagepath.disk,
                         distance.pixel = 117.9034,
                         known.distance = 1, low.size = 0.01)
chlor.leaf.area <- run.ij(path.imagej = ij.path,
                          set.directory = imagepath.chlLeaf,
                          distance.pixel = 117.9034,
                          known.distance = 1, low.size = 0.01)

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
  dplyr::select(id, abs.649, abs.665, chl.biomass = chlor.biomass, disk.area, chl.leaf.area) %>%
  mutate(chlA.ugml = (12.19 * abs.665) - (3.56 * abs.649),
         chlB.ugml = (21.99 * abs.649) - (5.32 * abs.665),
         chlA.ugml = ifelse(chlA.ugml < 0, 0, chlA.ugml),
         chlB.ugml = ifelse(chlB.ugml < 0, 0 , chlB.ugml),
         chlA.gml = chlA.ugml / 1000000,
         chlB.gml = chlB.ugml / 1000000,
         chlA.g = chlA.gml * 10, # extracted in 10mL DMSO
         chlB.g = chlB.gml * 10, # extracted in 10mL DMSO
         chlA.gm2 = chlA.g / (disk.area / 10000),
         chlB.gm2 = chlB.g / (disk.area / 10000),
         chlA.mmolm2 = chlA.gm2 / 893.51 * 1000,
         chlB.mmolm2 = chlB.gm2 / 907.47 * 1000,
         chl.marea = chl.biomass / (chl.leaf.area / 10000),
         chlA.mmolg = chlA.mmolm2 / chl.marea,
         chlB.mmolg = chlB.mmolm2 / chl.marea,
         chl.mmolg = chlA.mmolg + chlB.mmolg,
         chl.mmolm2 = chlA.mmolm2 + chlB.mmolm2,
         chlA.chlB = chlA.g / chlB.g) %>%
  dplyr::select(id, chlA.ugml:chlA.chlB)

###############################################################################
## Focal leaf area, leaf nitrogen content
###############################################################################
## Load CN data files
cn.files <- list.files(path = "../data_sheets/costech_results/", 
                          pattern = "\\.csv$", recursive = TRUE,
                          full.names = TRUE)
cn.files <- setNames(cn.files, cn.files)

cn.data <- lapply(cn.files, read.csv) %>% reshape::merge_all() %>%
  filter(type == "unknown") %>%
  mutate(nmass = as.numeric(n.weight.percent) / 100,
         cmass = as.numeric(c.weight.percent) / 100) %>%
  dplyr::select(id, nmass, cmass) %>%
  separate(id, into = c("co2", "inoc", "n.trt", "rep", "organ"), sep = "_") %>%
  unite("id", co2:rep, sep = "_") %>%
  pivot_wider(names_from = organ, values_from = nmass:cmass, names_sep = ".")

## Calculate leaf disk area
imagepath.focal <- "/Users/eaperkowski/git/2022_NxCO2xI/leaf_area/focal_scans/"

focal.area <- run.ij(path.imagej = ij.path,
                     set.directory = imagepath.focal,
                     distance.pixel = 117.9034,
                     known.distance = 1, low.size = 0.05,
                     set.memory = 30)

names(focal.area)[1:2] <- c("id", "focal.area")
focal.area$id <- gsub("_focal", "", focal.area$id)

###############################################################################
## %Ndfa
###############################################################################
head(isotopes)

df.15n <- isotopes %>%
  separate(id, into = c("co2", "inoc", "n.trt", "rep"), remove = FALSE) %>%
  unite("trt.combs", co2:n.trt, remove = FALSE) %>%
  full_join(biomass_area, by = "id") %>%
  dplyr::select(id:rep, leaf.d15n, nodule.biomass)

# Isolate B value (individuals totally reliant on Nfixation)
b.15n <- df.15n %>%
  filter(trt.combs == "e_y_0" | trt.combs == "a_y_0") %>%
  group_by(co2) %>%
  summarize(B = mean(leaf.d15n, na.rm = TRUE))

# Calculate ref.15N value for each CO2 x N combination for uninoculated pots
ref.15n <- df.15n %>%
  filter(inoc != "y" & nodule.biomass < 0.05) %>%
  filter(id != "a_n_0_111") %>%
  group_by(co2, n.trt) %>%
  summarize(ref.15n = mean(leaf.d15n, na.rm = TRUE))

d15n.merged <- df.15n %>%
  full_join(b.15n) %>%
  full_join(ref.15n)

ndfa <- d15n.merged %>%
  mutate(ndfa = calc_ndfa(ref.15n = ref.15n, 
                          sample.15n = leaf.d15n, B = B),
         ndfa = ifelse(ndfa > 100, 
                       100, 
                       ifelse(ndfa < 0, 0, ndfa))) %>%
  dplyr::select(id, leaf.d15n, ndfa)


###############################################################################
## Calculate beta
###############################################################################
head(isotopes)

d13c_air$co2_cat <- as.double(d13c_air$co2_cat)


beta <- isotopes %>%
  dplyr::select(id, leaf.d13c) %>%
  separate(id, into = c("co2", "inoc", "n.trt", "rep"), remove = FALSE) %>%
  mutate(temp = ifelse(co2 == "e", 21.5, 21.3),
         rh = ifelse(co2 == "e", 51.6, 50.3),
         co2_cat = ifelse(co2 == "e", 1000, 420),
         co2_num = ifelse(co2 == "e", 989, 439),
         vpd = RHtoVPD(RH = rh, TdegC = temp)) %>%
  full_join(d13c_air) %>%
  select(-co2.ppm) %>%
  mutate(chi = calc_chi_c3(leaf.d13c = leaf.d13c, air = d13c.air)[[2]]) %>%
  mutate(beta = calc_beta(chi = chi, temp = temp, vpd = vpd * 1000, 
                          ca = co2_num, z = 976)$beta) %>%
  dplyr::select(id, chi, beta)

###############################################################################
## Compile data files into single file for analyses/figs
###############################################################################
compile_df <- id %>%
  full_join(focal.area) %>% 
  full_join(biomass_area) %>%
  full_join(cn.data) %>%
  full_join(chlorophyll) %>%
  full_join(chl.leaf.area) %>%
  full_join(photo) %>%
  full_join(photo.wk6) %>%
  full_join(isotopes) %>%
  full_join(ndfa) %>%
  full_join(beta) %>%
  separate(id, into = c("co2", "inoc", "n.trt", "rep"), remove = FALSE) %>%
  mutate(rep = str_pad(rep, width = 3, side = "left", pad = "0"),
         
         ## Leaf N content
         marea = focal.biomass / (focal.area / 10000),
         marea.chl = chlor.biomass / (chl.leaf.area / 10000),
         narea = nmass.focal * marea,
         narea.chl = nmass.focal * marea.chl,
         
         ## Proportion of N calculations
         p.rubisco = p_rubisco(vcmax25, narea),
         p.bioe = p_bioenergetics(jmax25, narea),
         p.lightharv = p_lightharvesting(chl.mmolg, nmass.focal), ## swap with nmass.chl once data are in
         p.photo = p.rubisco + p.bioe + p.lightharv,
         p.structure = p_structure(lma = marea, narea = narea),
         
         ## Nitrogen-water use tradeoffs
         pnue = anet / (narea / 14),
         iwue = anet / gsw,
         narea.chi = narea / chi,
         vcmax.chi = vcmax25 / chi,
         
         ## stomatal limitation
         stomlim = stomatal_limitation(A_net = anet, Vcmax = vcmax25, leaf.temp = tleaf,
                                       Rd.meas = TRUE, Rd = rd25, temp = "C")[[5]],
         
         ## Tissue C and N biomass. Note that chlorophyll biomass is multiplied
         ## by nmass.focal due to high correlation between both
         leaf.totaln = (nmass.tl * leaf.biomass) + 
           (nmass.focal * focal.biomass) + (nmass.focal * chlor.biomass),
         stem.totaln = nmass.ts * stem.biomass,
         root.totaln = nmass.tr * root.biomass,
         root.totalc = cmass.tr * root.biomass,
         nod.totaln = ifelse(is.na(nmass.nod),
                             0, nmass.nod * nodule.biomass),
         nod.totalc = ifelse(is.na(cmass.nod),
                             0, cmass.nod * nodule.biomass),
         
         ## Ncost calcs
         wpn = leaf.totaln + stem.totaln + root.totaln + nod.totaln,
         cbg = root.totalc + nod.totalc,
         ncost = cbg / wpn,
         
         ## Whole plant growth
         tla.full = tla + focal.area + chl.leaf.area + disk.area,
         tla.full = ifelse(id == "e_y_70_10",
                           tla + focal.area + focal.area, tla.full),
         tla = tla.full,
         nodule.biomass = ifelse(inoc == "n" & is.na(nodule.biomass),
                                 0, nodule.biomass)) %>%
  arrange(rep) %>%
  dplyr::select(-tla.full, -notes) %>%
  mutate(across(.cols = c(nmass.focal:chlB.ugml,
                          marea:ncost),
         round, digits = 4)) %>%
  mutate(co2 = ifelse(co2 == "a", "amb", "elv"),
         inoc = ifelse(inoc == "n", "no.inoc", "inoc")) %>%
  as.data.frame()

write_csv(compile_df, "../data_sheets/NxCO2xI_compiled_datasheet.csv")




