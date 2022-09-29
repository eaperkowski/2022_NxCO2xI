# Compiling script for CO2 response curves and dark respiration values
# Note that all paths in the first part of this script assume working directory 
# root is the path to this script in the "2022_NxCO2_growthchamber" git
# repository. Path assigned via Apple operating system, may differ  on Windows 
# operating systems.

###############################################################################
## Load readLicorData package
###############################################################################
library(readLicorData)
library(dplyr)

###############################################################################
## Set working directory to Brad flash drive ("BP") and Dinah licor folder
###############################################################################
setwd("")

###############################################################################
## Load pre heatwave measurements. Adding column to indicate whether
## curves were done before or after heat stress event
###############################################################################
# NOTE: All files from first week of measurements were concatenated,
#       ("appended" on the machine). Object (e.g., "preHeat_alb") should
#       therefore contain data from both Sept. 7 and 8

# Albert
preHeat_alb <- licorData(location = "Albert/dinah_potato/2022-09-07_dinah_potato_n") %>%
  mutate(meas.type = "pre_heatwave") %>%
  slice(-(952:1016)) ## Remove appended file header from Sept. 8
# write.csv(preHeat_alb, "cleaned_files/Dinah_potato_preHeat_alb.csv")

# Stan
preHeat_stan <- licorData(location = "Stan/dinah_potato/2022-09-07_dinah_potato_n") %>%
  mutate(meas.type = "pre_heatwave") %>%
  slice(-(858:922)) ## Remove appended file header from Sept. 8
# write.csv(preHeat_stan, "cleaned_files/Dinah_potato_preHeat_stan.csv")

# Gibson
preHeat_gib <- licorData(location = "Gibson/2022-09-07-Dinah-potato/2022-09-07-Dinah Nitrogen") %>%
  mutate(meas.type = "pre_heatwave") %>%
  slice(-(714:778)) ## Remove appended file header from Sept. 8
# NOTE: Gibson experienced flow leaks on day 2 and stopped taking measurements 
#       early, which is why there are fewer observations in file
# write.csv(preHeat_gib, "cleaned_files/Dinah_potato_preHeat_gib.csv")

# Ozzie
preHeat_ozz <- licorData(location = "Ozzie/2022-09-07-Dinah-Potato/2022-09-07-0945_potato_N") %>%
  mutate(meas.type = "pre_heatwave") %>%
  slice(-(1233:1297)) ## Remove appended file header from Sept. 8
# write.csv(preHeat_ozz, "cleaned_files/Dinah_potato_preHeat_ozz.csv")

###############################################################################
## Load post heatwave measurements. Adding column to indicate whether
## curves were done before or after heat stress event
###############################################################################

# Albert
postHeat_alb_d1 <- licorData(location = "Albert/2022-09-15-dinah-potato-heatday1/2022-09-15-0854_dinah-potato-heatd1") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_alb_d1, "cleaned_files/Dinah_potato_postHeat_alb_d1.csv")

postHeat_alb_d2 <- licorData(location = "Albert/dinah_potato/2022-09-16_temperature_day2") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_alb_d2, "cleaned_files/Dinah_potato_postHeat_alb_d2.csv")

# Stan
postHeat_stan_d1 <- licorData(location = "Stan/2022-09-15-dinah-potato-heatday1/2022-09-15-0855_dinah-potato-heatday1") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_stan_d1, "cleaned_files/Dinah_potato_postHeat_stan_d1.csv")

postHeat_stan_d2 <- licorData(location = "Stan/dinah_potato/2022-09-16-0852_dinah_potato_heatday2") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_stan_d2, "cleaned_files/Dinah_potato_postHeat_stan_d2.csv")

# Gibson
postHeat_gib_d1 <- licorData(location = "Gibson/2022-09-07-Dinah-potato/2022-09-15-dinah-potato-heatday1/2022-09-15-0939_dinah-potato-heatday1") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_gib_d1, "cleaned_files/Dinah_potato_postHeat_gib_d1.csv")

postHeat_gib_d2 <- licorData(location = "Gibson/2022-09-07-Dinah-potato/2022-09-16-temperature_day2") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_gib_d2, "cleaned_files/Dinah_potato_postHeat_gib_d2.csv")

# Ozzie
postHeat_ozz_d1 <- licorData(location = "Ozzie/2022-09-15-Dinah-potato-heatday1/2022-09-15-0852_dinah_heat_d1") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_ozz_d1, "cleaned_files/Dinah_potato_postHeat_ozz_d1.csv")

postHeat_ozz_d2 <- licorData(location = "Ozzie/2022-09-07-Dinah-Potato/2022-09-16-0851_dinah_potato_heatday2") %>%
  mutate(meas.type = "post_heatwave")
# write.csv(postHeat_ozz_d2, "cleaned_files/Dinah_potato_postHeat_ozz_d2.csv")


###############################################################################
## Merge files into central file. Useful for 'fitacis' when fitting multiple
## curves
###############################################################################
# NOTE: Using list.files notation to avoid common merge conflict with 
# readLicorData package. Cols seem to be assigned different classes when
# cleaned through 'licorData', which makes merging files difficult/unnecessarily
# time consuming. Reloading files into list of data frames, them merging through
# reshape::merge_all() seems to do the trick.

# List files
file.list <- list.files("cleaned_files",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list <- setNames(file.list, stringr::str_extract(basename(file.list), 
                                                      '.*(?=\\.csv)'))

# Merge list of data frames, arrange by marchine, measurement type, id, and time elapsed
merged_curves <- lapply(file.list, read.csv) %>%
  reshape::merge_all() %>%
  arrange(machine, meas.type, id, elapsed)
#write.csv(merged_curves, "./Dinah_potato_curves_fullyMerged.csv")

## End of data cleaning, ready for curve fitting ##