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
## Load co2 response curve data files for week 6. Adding column to indicate whether
## curves were done on week 6 or week 7 of growth
###############################################################################

# Albert
eco2_co2Resp_w6d1_alb <- licorData(location = "../licor_raw/week6/eco2_co2Resp_W6D1_alb") %>%
  mutate(week = 6)
# write.csv(eco2_co2Resp_w6d1_alb, "../licor_cleaned/co2_resp/eco2_co2Resp_w6d1_alb.csv", row.names = FALSE)

eco2_co2Resp_w6d2_alb <- licorData(location = "../licor_raw/week6/eco2_co2Resp_W6D2_alb") %>%
  mutate(week = 6)
# write.csv(eco2_co2Resp_w6d2_alb, "../licor_cleaned/co2_resp/eco2_co2Resp_w6d2_alb.csv", row.names = FALSE)

aco2_co2Resp_w6d1_alb <- licorData(location = "../licor_raw/week6/aco2_co2Resp_W6D1_alb") %>%
  mutate(week = 6)
# write.csv(aco2_co2Resp_w6d1_alb, "../licor_cleaned/co2_resp/aco2_co2Resp_w6d1_alb.csv", row.names = FALSE)

aco2_co2Resp_w6d2_alb <- licorData(location = "../licor_raw/week6/aco2_co2Resp_W6D2_alb") %>%
  mutate(week = 6)
# write.csv(aco2_co2Resp_w6d2_alb, "../licor_cleaned/co2_resp/aco2_co2Resp_w6d2_alb.csv", row.names = FALSE)

# Stan
eco2_co2Resp_w6d1_stan <- licorData(location = "../licor_raw/week6/eco2_co2Resp_W6D1_stan") %>%
  mutate(week = 6)
# write.csv(eco2_co2Resp_w6d1_stan, "../licor_cleaned/co2_resp/eco2_co2Resp_w6d1_stan.csv", row.names = FALSE)

eco2_co2Resp_w6d2_stan <- licorData(location = "../licor_raw/week6/eco2_co2Resp_W6D2_stan") %>%
  mutate(week = 6)
# write.csv(eco2_co2Resp_w6d2_stan, "../licor_cleaned/co2_resp/eco2_co2Resp_w6d2_stan.csv", row.names = FALSE)

aco2_co2Resp_w6d1_stan <- licorData(location = "../licor_raw/week6/aco2_co2Resp_W6D1_stan") %>%
  mutate(week = 6)
# write.csv(aco2_co2Resp_w6d1_stan, "../licor_cleaned/co2_resp/aco2_co2Resp_w6d1_stan.csv", row.names = FALSE)

aco2_co2Resp_w6d2_stan <- licorData(location = "../licor_raw/week6/aco2_co2Resp_W6D2_stan") %>%
  mutate(week = 6)
# write.csv(aco2_co2Resp_w6d2_stan, "../licor_cleaned/co2_resp/aco2_co2Resp_w6d2_stan.csv", row.names = FALSE)

# Ozzie
eco2_co2Resp_w6d1_ozz <- licorData(location = "../licor_raw/week6/eco2_co2Resp_W6D1_ozz") %>%
  mutate(week = 6)
# write.csv(eco2_co2Resp_w6d1_ozz, "../licor_cleaned/co2_resp/eco2_co2Resp_w6d1_ozz.csv", row.names = FALSE)

eco2_co2Resp_w6d2_ozz <- licorData(location = "../licor_raw/week6/eco2_co2Resp_W6D2_ozz") %>%
  mutate(week = 6)
# write.csv(eco2_co2Resp_w6d2_ozz, "../licor_cleaned/co2_resp/eco2_co2Resp_w6d2_ozz.csv", row.names = FALSE)

aco2_co2Resp_w6d1_ozz <- licorData(location = "../licor_raw/week6/aco2_co2Resp_W6D1_ozz") %>%
  mutate(week = 6)
# write.csv(aco2_co2Resp_w6d1_ozz, "../licor_cleaned/co2_resp/aco2_co2Resp_w6d1_ozz.csv", row.names = FALSE)

aco2_co2Resp_w6d2_ozz <- licorData(location = "../licor_raw/week6/aco2_co2Resp_W6D2_ozz") %>%
  mutate(week = 6)
# write.csv(aco2_co2Resp_w6d2_ozz, "../licor_cleaned/co2_resp/aco2_co2Resp_w6d2_ozz.csv", row.names = FALSE)

###############################################################################
## Load co2 response curve data files for week 7. Adding column to indicate 
## whether curves were done on week 6 or week 7 of growth
###############################################################################

# Albert
eco2_co2Resp_w7d1_alb <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D1_alb") %>%
  mutate(week = 7)
# write.csv(eco2_co2Resp_w7d1_alb, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d1_alb.csv", row.names = FALSE)

eco2_co2Resp_w7d2_alb <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D2_alb") %>%
  mutate(week = 7)
# write.csv(eco2_co2Resp_w7d2_alb, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d2_alb.csv", row.names = FALSE)

aco2_co2Resp_w7d1_alb <- licorData(location = "../licor_raw/week7/aco2_co2Resp_W7D1_alb") %>%
  mutate(week = 7)
# write.csv(aco2_co2Resp_w7d1_alb, "../licor_cleaned/co2_resp/aco2_co2Resp_w7d1_alb.csv", row.names = FALSE)

aco2_co2Resp_w7d2_alb <- licorData(location = "../licor_raw/week7/aco2_co2Resp_W7D2_alb") %>%
  mutate(week = 7)
# write.csv(aco2_co2Resp_w7d2_alb, "../licor_cleaned/co2_resp/aco2_co2Resp_w7d2_alb.csv", row.names = FALSE)

# Stan
eco2_co2Resp_w7d1_stan <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D1_stan") %>%
  mutate(week = 7)
# write.csv(eco2_co2Resp_w7d1_stan, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d1_stan.csv", row.names = FALSE)

eco2_co2Resp_w7d2_stan <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D2_stan") %>%
  mutate(week = 7)
# write.csv(eco2_co2Resp_w7d2_stan, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d2_stan.csv", row.names = FALSE)

aco2_co2Resp_w7d1_stan <- licorData(location = "../licor_raw/week7/aco2_co2Resp_W7D1_stan") %>%
  mutate(week = 7)
# write.csv(aco2_co2Resp_w7d1_stan, "../licor_cleaned/co2_resp/aco2_co2Resp_w7d1_stan.csv", row.names = FALSE)

aco2_co2Resp_w7d2_stan <- licorData(location = "../licor_raw/week7/aco2_co2Resp_W7D2_stan") %>%
  mutate(week = 7)
# write.csv(aco2_co2Resp_w7d2_stan, "../licor_cleaned/co2_resp/aco2_co2Resp_w7d2_stan.csv", row.names = FALSE)

# Ozzie
eco2_co2Resp_w7d1_ozz <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D1_ozz") %>%
  mutate(week = 7)
# write.csv(eco2_co2Resp_w7d1_ozz, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d1_ozz.csv", row.names = FALSE)

eco2_co2Resp_w7d2_ozz <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D2_ozz") %>%
  mutate(week = 7)
# write.csv(eco2_co2Resp_w7d2_ozz, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d2_ozz.csv", row.names = FALSE)

aco2_co2Resp_w7d1_ozz <- licorData(location = "../licor_raw/week7/aco2_co2Resp_W7D1_ozz") %>%
  mutate(week = 7)
# write.csv(aco2_co2Resp_w7d1_ozz, "../licor_cleaned/co2_resp/aco2_co2Resp_w7d1_ozz.csv", row.names = FALSE)

aco2_co2Resp_w7d2_ozz <- licorData(location = "../licor_raw/week7/aco2_co2Resp_W7D2_ozz") %>%
  mutate(week = 7)
# write.csv(aco2_co2Resp_w7d2_ozz, "../licor_cleaned/co2_resp/aco2_co2Resp_w7d2_ozz.csv", row.names = FALSE)

###############################################################################
## Load dark respiration data files for week 6. Adding column to indicate 
## whether curves were done on week 6 or week 7 of growth
###############################################################################

# Albert
eco2_rd_w6d1_alb <- licorData(location = "../licor_raw/week6/eco2_rd_W6D1_alb") %>%
  mutate(week = 6)
write.csv(eco2_rd_w6d1_alb, "../licor_cleaned/rd/eco2_rd_w6d1_alb.csv", row.names = FALSE)

eco2_rd_w6d2_alb <- licorData(location = "../licor_raw/week6/eco2_rd_W6D2_alb") %>%
  mutate(week = 6)
write.csv(eco2_rd_w6d2_alb, "../licor_cleaned/rd/eco2_rd_w6d2_alb.csv", row.names = FALSE)

aco2_rd_w6d1_alb <- licorData(location = "../licor_raw/week6/aco2_rd_W6D1_alb") %>%
  mutate(week = 6)
write.csv(aco2_rd_w6d1_alb, "../licor_cleaned/rd/aco2_rd_w6d1_alb.csv", row.names = FALSE)

aco2_rd_w6d2_alb <- licorData(location = "../licor_raw/week6/aco2_rd_W6D2_alb") %>%
  mutate(week = 6)
write.csv(aco2_rd_w6d2_alb, "../licor_cleaned/rd/aco2_rd_w6d2_alb.csv", row.names = FALSE)

# Stan
eco2_rd_w6d1_stan <- licorData(location = "../licor_raw/week6/eco2_rd_W6D1_stan") %>%
  mutate(week = 6)
write.csv(eco2_rd_w6d1_stan, "../licor_cleaned/rd/eco2_rd_w6d1_stan.csv", row.names = FALSE)

eco2_rd_w6d2_stan <- licorData(location = "../licor_raw/week6/eco2_rd_W6D2_stan") %>%
  mutate(week = 6)
write.csv(eco2_rd_w6d2_stan, "../licor_cleaned/rd/eco2_rd_w6d2_stan.csv", row.names = FALSE)

aco2_rd_w6d1_stan <- licorData(location = "../licor_raw/week6/aco2_rd_W6D1_stan") %>%
  mutate(week = 6)
write.csv(aco2_rd_w6d1_stan, "../licor_cleaned/rd/aco2_rd_w6d1_stan.csv", row.names = FALSE)

aco2_rd_w6d2_stan <- licorData(location = "../licor_raw/week6/aco2_rd_W6D2_stan") %>%
  mutate(week = 6)
write.csv(aco2_rd_w6d2_stan, "../licor_cleaned/rd/aco2_rd_w6d2_stan.csv", row.names = FALSE)

# Ozzie
eco2_rd_w6d1_ozz <- licorData(location = "../licor_raw/week6/eco2_rd_W6D1_ozz") %>%
  mutate(week = 6)
write.csv(eco2_rd_w6d1_ozz, "../licor_cleaned/rd/eco2_rd_w6d1_ozz.csv", row.names = FALSE)

eco2_rd_w6d2_ozz <- licorData(location = "../licor_raw/week6/eco2_rd_W6D2_ozz") %>%
  mutate(week = 6)
write.csv(eco2_rd_w6d2_ozz, "../licor_cleaned/rd/eco2_rd_w6d2_ozz.csv", row.names = FALSE)

aco2_rd_w6d1_ozz <- licorData(location = "../licor_raw/week6/aco2_rd_W6D1_ozz") %>%
  mutate(week = 6)
write.csv(aco2_rd_w6d1_ozz, "../licor_cleaned/rd/aco2_rd_w6d1_ozz.csv", row.names = FALSE)

aco2_rd_w6d2_ozz <- licorData(location = "../licor_raw/week6/aco2_rd_W6D2_ozzie") %>%
  mutate(week = 6)
write.csv(aco2_rd_w6d2_ozz, "../licor_cleaned/rd/aco2_rd_w6d2_ozz.csv", row.names = FALSE)

###############################################################################
## Load dark respiration data files for week 7. Adding column to indicate 
## whether curves were done on week 6 or week 7 of growth
###############################################################################

# Albert
eco2_rd_w7d1_alb <- licorData(location = "../licor_raw/week7/eco2_rd_W7D1_alb") %>%
  mutate(week = 7)
write.csv(eco2_rd_w7d1_alb, "../licor_cleaned/rd/eco2_rd_w7d1_alb.csv", row.names = FALSE)

eco2_rd_w7d2_alb <- licorData(location = "../licor_raw/week7/eco2_rd_W7D2_alb") %>%
  mutate(week = 7)
write.csv(eco2_rd_w7d2_alb, "../licor_cleaned/rd/eco2_rd_w7d2_alb.csv", row.names = FALSE)

aco2_rd_w7d1_alb <- licorData(location = "../licor_raw/week7/aco2_rd_W7D1_alb") %>%
  mutate(week = 7)
write.csv(aco2_rd_w7d1_alb, "../licor_cleaned/rd/aco2_rd_w7d1_alb.csv", row.names = FALSE)

aco2_rd_w7d2_alb <- licorData(location = "../licor_raw/week7/aco2_rd_W7D2_alb") %>%
  mutate(week = 7)
write.csv(aco2_rd_w7d2_alb, "../licor_cleaned/rd/aco2_rd_w7d2_alb.csv", row.names = FALSE)

# Stan
eco2_rd_w7d1_stan <- licorData(location = "../licor_raw/week7/eco2_rd_W7D1_stan") %>%
  mutate(week = 7)
write.csv(eco2_rd_w7d1_stan, "../licor_cleaned/rd/eco2_rd_w7d1_stan.csv", row.names = FALSE)

eco2_rd_w7d2_stan <- licorData(location = "../licor_raw/week7/eco2_rd_W7D2_stan") %>%
  mutate(week = 7)
write.csv(eco2_rd_w7d2_stan, "../licor_cleaned/rd/eco2_rd_w7d2_stan.csv", row.names = FALSE)

aco2_rd_w7d1_stan <- licorData(location = "../licor_raw/week7/aco2_rd_W7D1_stan") %>%
  mutate(week = 7)
write.csv(aco2_rd_w7d1_stan, "../licor_cleaned/rd/aco2_rd_w7d1_stan.csv", row.names = FALSE)

aco2_rd_w7d2_stan <- licorData(location = "../licor_raw/week7/aco2_rd_W7D2_stan") %>%
  mutate(week = 7)
write.csv(aco2_rd_w7d2_stan, "../licor_cleaned/rd/aco2_rd_w7d2_stan.csv", row.names = FALSE)

# Ozzie
eco2_rd_w7d1_ozz <- licorData(location = "../licor_raw/week7/eco2_rd_W7D1_ozz") %>%
  mutate(week = 7)
write.csv(eco2_rd_w7d1_ozz, "../licor_cleaned/rd/eco2_rd_w7d1_ozz.csv", row.names = FALSE)

eco2_rd_w7d2_ozz <- licorData(location = "../licor_raw/week7/eco2_rd_W7D2_ozz") %>%
  mutate(week = 7)
write.csv(eco2_rd_w7d2_ozz, "../licor_cleaned/rd/eco2_rd_w7d2_ozz.csv", row.names = FALSE)

aco2_rd_w7d1_ozz <- licorData(location = "../licor_raw/week7/aco2_rd_W7D1_ozz") %>%
  mutate(week = 7)
write.csv(aco2_rd_w7d1_ozz, "../licor_cleaned/rd/aco2_rd_w7d1_ozz.csv", row.names = FALSE)

aco2_rd_w7d2_ozz <- licorData(location = "../licor_raw/week7/aco2_rd_W7D2_ozz") %>%
  mutate(week = 7)
write.csv(aco2_rd_w7d2_ozz, "../licor_cleaned/rd/aco2_rd_w7d2_ozz.csv", row.names = FALSE)

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
file.list <- list.files("../licor_cleaned/co2_resp",
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