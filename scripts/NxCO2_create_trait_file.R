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
curve.fits <- read.csv("../data_sheets/NxCO2_curve_results.csv")
curve.fits[98,1] <- "a_y_280_98"

###############################################################################
## Load propN functions
###############################################################################
source("../../r_functions/propN_funcs.R")

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
  select(id, abs.649, abs.665, chl.biomass = chlor.biomass, disk.area, chl.leaf.area) %>%
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
  select(id, chlA.ugml:chlA.chlB)

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
  select(id, nmass, cmass) %>%
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

## Compile data file into single file (to be used for stats/figs)
compile_df <- focal.area %>% 
  full_join(biomass_area) %>%
  full_join(cn.data) %>%
  full_join(chlorophyll) %>%
  full_join(chl.leaf.area) %>%
  full_join(curve.fits) %>%
  separate(id, into = c("co2", "inoc", "n.trt", "rep"), remove = FALSE) %>%
  mutate(rep = str_pad(rep, width = 3, side = "left", pad = "0"),
         
         ## Leaf N content
         marea = focal.biomass / (focal.area / 10000),
         marea.chl = chlor.biomass / (chl.leaf.area / 10000),
         narea = nmass.focal * marea,
         #narea.chl = nmass.chl * marea.chl,
         
         
         ## Proportion of N calculations
         p.rubisco = p_rubisco(vcmax25, narea),
         p.bioe = p_bioenergetics(jmax25, narea),
         p.lightharv = p_lightharvesting(chl.mmolg, nmass.focal), ## swap with nmass.chl once data are in
         p.photo = p.rubisco + p.bioe + p.lightharv,
         #p.structure = p_structure(lma = marea, narea = narea, useEq = FALSE),
         
         ## Tissue C and N biomasses
         leaf.totaln = nmass.tl * leaf.biomass, #+ nmass.nod * chor.biomass,
         stem.totaln = nmass.ts * stem.biomass,
         root.totaln = nmass.tr * root.biomass,
         root.totalc = cmass.tr * root.biomass,
         nod.totaln = ifelse(nodule.biomass == 0 | is.na(nmass.nod),
                             0, nmass.nod * nodule.biomass),
         nod.totalc = ifelse(nodule.biomass == 0 | is.na(cmass.nod),
                             0, cmass.nod * nodule.biomass),
         
         ## Ncost calcs
         wpn = leaf.totaln + stem.totaln + root.totaln + nod.totaln,
         cbg = root.totalc + nod.totalc,
         ncost = cbg / wpn,
         
         ## Whole plant growth
         tla = tla + focal.area + chl.leaf.area + disk.area,
         nodule.biomass = ifelse(inoc == "n" & is.na(nodule.biomass),
                                 0, nodule.biomass)) %>%
  arrange(rep)

write.csv(compile_df,
          "../data_sheets/NxCO2xI_compiled_datasheet.csv", row.names = FALSE)


ggplot(data = subset(compile_df, marea.chl < 100), 
       aes(x = marea, y = marea.chl)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', size = 1) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  scale_x_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  scale_y_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  labs(x = expression("M"["area"]*" (g m"^"-2"*")"),
       y = expression("M"["area_chl"]*" (g m"^"-2"*")")) +
  theme_bw(base_size = 18)


ggplot(data = subset(compile_df, co2 == "e"), 
       aes(x = as.numeric(n.trt), y = ncost, fill = inoc)) +
  geom_point(shape = 21, size = 4, alpha = 0.75) +
  scale_fill_discrete(labels = c("no", "yes")) +
  geom_smooth(method = 'lm', formula = y ~ poly(x, 2)) +
  labs(x = "Soil nitrogen fertilization (ppm)",
       y = expression("Carbon cost to acquire nitrogen (gC gN"^-1*")"),
       fill = "Inoculation status") +
  theme_bw(base_size = 18)




ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = narea, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = vcmax25, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = jmax25, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = p.photo, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = p.bioe, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = p.lightharv, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = marea, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = nmass.focal, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = nodule.biomass, fill = co2)) +
  geom_point(shape = 21, size = 4, alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = nodule.biomass / root.biomass, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = total.biomass, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)

ggplot(data = compile_df, aes(x = as.numeric(n.trt), 
                              y = tla, fill = co2)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  facet_grid(~inoc)
