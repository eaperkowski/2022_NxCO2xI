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

eco2_co2Resp_w7d3_ozz <- licorData(location = "../licor_raw/week7/eco2_co2Resp_W7D3_ozz") %>%
  mutate(week = 7)
#write.csv(eco2_co2Resp_w7d3_ozz, "../licor_cleaned/co2_resp/eco2_co2Resp_w7d3_ozz.csv", row.names = FALSE)


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
# write.csv(eco2_rd_w6d1_alb, "../licor_cleaned/rd/eco2_rd_w6d1_alb.csv", row.names = FALSE)

eco2_rd_w6d2_alb <- licorData(location = "../licor_raw/week6/eco2_rd_W6D2_alb") %>%
  mutate(week = 6)
# write.csv(eco2_rd_w6d2_alb, "../licor_cleaned/rd/eco2_rd_w6d2_alb.csv", row.names = FALSE)

aco2_rd_w6d1_alb <- licorData(location = "../licor_raw/week6/aco2_rd_W6D1_alb") %>%
  mutate(week = 6)
# write.csv(aco2_rd_w6d1_alb, "../licor_cleaned/rd/aco2_rd_w6d1_alb.csv", row.names = FALSE)

aco2_rd_w6d2_alb <- licorData(location = "../licor_raw/week6/aco2_rd_W6D2_alb") %>%
  mutate(week = 6)
# write.csv(aco2_rd_w6d2_alb, "../licor_cleaned/rd/aco2_rd_w6d2_alb.csv", row.names = FALSE)

# Stan
eco2_rd_w6d1_stan <- licorData(location = "../licor_raw/week6/eco2_rd_W6D1_stan") %>%
  mutate(week = 6)
# write.csv(eco2_rd_w6d1_stan, "../licor_cleaned/rd/eco2_rd_w6d1_stan.csv", row.names = FALSE)

eco2_rd_w6d2_stan <- licorData(location = "../licor_raw/week6/eco2_rd_W6D2_stan") %>%
  mutate(week = 6)
# write.csv(eco2_rd_w6d2_stan, "../licor_cleaned/rd/eco2_rd_w6d2_stan.csv", row.names = FALSE)

aco2_rd_w6d1_stan <- licorData(location = "../licor_raw/week6/aco2_rd_W6D1_stan") %>%
  mutate(week = 6)
# write.csv(aco2_rd_w6d1_stan, "../licor_cleaned/rd/aco2_rd_w6d1_stan.csv", row.names = FALSE)

aco2_rd_w6d2_stan <- licorData(location = "../licor_raw/week6/aco2_rd_W6D2_stan") %>%
  mutate(week = 6)
# write.csv(aco2_rd_w6d2_stan, "../licor_cleaned/rd/aco2_rd_w6d2_stan.csv", row.names = FALSE)

# Ozzie
eco2_rd_w6d1_ozz <- licorData(location = "../licor_raw/week6/eco2_rd_W6D1_ozz") %>%
  mutate(week = 6)
# write.csv(eco2_rd_w6d1_ozz, "../licor_cleaned/rd/eco2_rd_w6d1_ozz.csv", row.names = FALSE)

eco2_rd_w6d2_ozz <- licorData(location = "../licor_raw/week6/eco2_rd_W6D2_ozz") %>%
  mutate(week = 6)
# write.csv(eco2_rd_w6d2_ozz, "../licor_cleaned/rd/eco2_rd_w6d2_ozz.csv", row.names = FALSE)

aco2_rd_w6d1_ozz <- licorData(location = "../licor_raw/week6/aco2_rd_W6D1_ozz") %>%
  mutate(week = 6)
# write.csv(aco2_rd_w6d1_ozz, "../licor_cleaned/rd/aco2_rd_w6d1_ozz.csv", row.names = FALSE)

aco2_rd_w6d2_ozz <- licorData(location = "../licor_raw/week6/aco2_rd_W6D2_ozzie") %>%
  mutate(week = 6)
# write.csv(aco2_rd_w6d2_ozz, "../licor_cleaned/rd/aco2_rd_w6d2_ozz.csv", row.names = FALSE)

###############################################################################
## Load dark respiration data files for week 7. Adding column to indicate 
## whether curves were done on week 6 or week 7 of growth
###############################################################################

# Albert
eco2_rd_w7d1_alb <- licorData(location = "../licor_raw/week7/eco2_rd_W7D1_alb") %>%
  mutate(week = 7)
# write.csv(eco2_rd_w7d1_alb, "../licor_cleaned/rd/eco2_rd_w7d1_alb.csv", row.names = FALSE)

eco2_rd_w7d2_alb <- licorData(location = "../licor_raw/week7/eco2_rd_W7D2_alb") %>%
  mutate(week = 7)
# write.csv(eco2_rd_w7d2_alb, "../licor_cleaned/rd/eco2_rd_w7d2_alb.csv", row.names = FALSE)

aco2_rd_w7d1_alb <- licorData(location = "../licor_raw/week7/aco2_rd_W7D1_alb") %>%
  mutate(week = 7)
# write.csv(aco2_rd_w7d1_alb, "../licor_cleaned/rd/aco2_rd_w7d1_alb.csv", row.names = FALSE)

aco2_rd_w7d2_alb <- licorData(location = "../licor_raw/week7/aco2_rd_W7D2_alb") %>%
  mutate(week = 7)
# write.csv(aco2_rd_w7d2_alb, "../licor_cleaned/rd/aco2_rd_w7d2_alb.csv", row.names = FALSE)

# Stan
eco2_rd_w7d1_stan <- licorData(location = "../licor_raw/week7/eco2_rd_W7D1_stan") %>%
  mutate(week = 7)
# write.csv(eco2_rd_w7d1_stan, "../licor_cleaned/rd/eco2_rd_w7d1_stan.csv", row.names = FALSE)

eco2_rd_w7d2_stan <- licorData(location = "../licor_raw/week7/eco2_rd_W7D2_stan") %>%
  mutate(week = 7)
# write.csv(eco2_rd_w7d2_stan, "../licor_cleaned/rd/eco2_rd_w7d2_stan.csv", row.names = FALSE)

aco2_rd_w7d1_stan <- licorData(location = "../licor_raw/week7/aco2_rd_W7D1_stan") %>%
  mutate(week = 7)
# write.csv(aco2_rd_w7d1_stan, "../licor_cleaned/rd/aco2_rd_w7d1_stan.csv", row.names = FALSE)

aco2_rd_w7d2_stan <- licorData(location = "../licor_raw/week7/aco2_rd_W7D2_stan") %>%
  mutate(week = 7)
# write.csv(aco2_rd_w7d2_stan, "../licor_cleaned/rd/aco2_rd_w7d2_stan.csv", row.names = FALSE)

# Ozzie
eco2_rd_w7d1_ozz <- licorData(location = "../licor_raw/week7/eco2_rd_W7D1_ozz") %>%
  mutate(week = 7)
# write.csv(eco2_rd_w7d1_ozz, "../licor_cleaned/rd/eco2_rd_w7d1_ozz.csv", row.names = FALSE)

eco2_rd_w7d2_ozz <- licorData(location = "../licor_raw/week7/eco2_rd_W7D2_ozz") %>%
  mutate(week = 7)
# write.csv(eco2_rd_w7d2_ozz, "../licor_cleaned/rd/eco2_rd_w7d2_ozz.csv", row.names = FALSE)

aco2_rd_w7d1_ozz <- licorData(location = "../licor_raw/week7/aco2_rd_W7D1_ozz") %>%
  mutate(week = 7)
# write.csv(aco2_rd_w7d1_ozz, "../licor_cleaned/rd/aco2_rd_w7d1_ozz.csv", row.names = FALSE)

aco2_rd_w7d2_ozz <- licorData(location = "../licor_raw/week7/aco2_rd_W7D2_ozz") %>%
  mutate(week = 7)
# write.csv(aco2_rd_w7d2_ozz, "../licor_cleaned/rd/aco2_rd_w7d2_ozz.csv", row.names = FALSE)

###############################################################################
## Merge co2 response curves into single file. Useful for 'fitacis' when 
## fitting multiple curves
###############################################################################
# NOTE: Using list.files notation to avoid common merge conflict with 
# readLicorData package. Cols seem to be assigned different classes when
# cleaned through 'licorData', which makes merging files difficult/unnecessarily
# time consuming. Reloading files into list of data frames, then merging through
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
  arrange(machine, week, id, elapsed)

co2_resp_wk6 <- merged_curves %>%
  filter(week == 6)
write.csv(co2_resp_wk6, "../data_sheets/NxCO2_co2_resp_wk6.csv", row.names = FALSE)

co2_resp_wk7 <- merged_curves %>%
  filter(week == 7)
write.csv(co2_resp_wk7, "../data_sheets/NxCO2_co2_resp_wk7.csv", row.names = FALSE)

###############################################################################
## Merge rd files into single file. Useful for 'fitacis' when fitting multiple
## curves
###############################################################################
# NOTE: Using list.files notation to avoid common merge conflict with 
# readLicorData package. Cols seem to be assigned different classes when
# cleaned through 'licorData', which makes merging files difficult/unnecessarily
# time consuming. Reloading files into list of data frames, then merging through
# reshape::merge_all() seems to do the trick.

## Load temp_standardize fxn
source("/Users/eaperkowski/git/r_functions/temp_standardize.R")


# List files
file.list.rd <- list.files("../licor_cleaned/rd",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list.rd <- setNames(file.list.rd, stringr::str_extract(basename(file.list.rd), 
                                                      '.*(?=\\.csv)'))

# Merge list of data frames, arrange by machine, measurement type, id, and time elapsed
rd <- lapply(file.list.rd, read.csv) %>%
  reshape::merge_all()

rd.wk6 <- rd %>%
  filter(week == 6) %>%
  group_by(id, week) %>%
  mutate(A = ifelse(A > 0, NA, A),
         rd = abs(A)) %>%
  summarize(rd = mean(rd, na.rm = TRUE),
            tLeaf = mean(Tleaf, na.rm = TRUE),
            rd25 = temp_standardize(rd,
                                    estimate.type = "Rd",
                                    standard.to = 25,
                                    tLeaf = tLeaf,
                                    tGrow = 22.5,
                                    pft = "C3H")) %>%
  arrange(id)
write.csv(rd.wk6, "../data_sheets/NxCO2_rd_wk6.csv", row.names = FALSE)

rd.wk7 <- rd %>%
  filter(week == 7) %>%
  group_by(id, week) %>%
  mutate(A = ifelse(A > 0, NA, A),
         rd = abs(A)) %>%
  summarize(rd = mean(rd, na.rm = TRUE),
            tLeaf = mean(Tleaf, na.rm = TRUE),
            rd25 = temp_standardize(rd,
                                    estimate.type = "Rd",
                                    standard.to = 25,
                                    tLeaf = tLeaf,
                                    tGrow = 22.5,
                                    pft = "C3H")) %>%
  filter(id != "a_y_280_133") # Remove first iteration of rd; co2 cylinder ran out
rd.wk7$id[rd.wk7$id == "a_y_280_133_b"] <- "a_y_280_133"

write.csv(rd.wk7, "../data_sheets/NxCO2_rd_wk7.csv", row.names = FALSE)

## End of data cleaning, ready for curve fitting ##