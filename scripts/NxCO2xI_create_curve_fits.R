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

# For creating unique ID list
biomass_area <- read.csv("../data_sheets/NxCO2_tla_biomass_data.csv")
ids <- data.frame(id = biomass_area$id)
ids <- separate(ids, id, sep = "_", into = c("co2", "inoc", "n.trt", "rep"),
                remove = FALSE)

###############################################################################
## Load custom fxns
###############################################################################
source("/Users/eaperkowski/git/r_functions/temp_standardize.R")
source("/Users/eaperkowski/git/r_functions/stomatal_limitation.R")

###############################################################################
## Prep data frame to fit A/Ci curves
###############################################################################
aci.prep <- co2.response.wk7 %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair,	Tleaf) %>%
  arrange(id) %>%
  left_join(rd.wk7, by = "id") %>%
  dplyr::select(-week, -tLeaf) %>%
  group_by(id) %>%
  mutate(keep.row = "yes") %>%
  data.frame()

aci.temps <- aci.prep %>%
  group_by(id) %>%
  summarize(Tleaf = mean(Tleaf, na.rm = TRUE))

###############################################################################
## Run A/Ci curves
###############################################################################

## Remove rows based on A/Ci fits, and also include all points measured
## at 0 ppm CO2. Workshop w/ Licor noted that 0ppm CO2 turns off mixing fan. 
#aci.merged$keep.row[aci.merged$A < -1.5] <- "no"
aci.prep$keep.row[c()] <- "no"

#####################################################################
# A/Ci curves
#####################################################################

#######################################
# Ambient CO2 inoculated
#######################################
e_y_0_1 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_0_1" &
                                 Ci < 850) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_1)
aci.coefs <- data.frame(id = "e_y_0_1", t(coef(e_y_0_1)))


e_y_0_2 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_0_2") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_2)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_0_2", t(coef(e_y_0_2)))


e_y_0_3 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_0_3") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_3)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_0_3", t(coef(e_y_0_3)))


e_y_0_4 <- aci.prep %>% filter(keep.row == "yes" & 
                                 id == "e_y_0_4") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_4)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_0_4", t(coef(e_y_0_4)))


e_y_35_5 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_5") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_5)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_35_5", t(coef(e_y_35_5)))


e_y_35_6 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_6") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_6)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_35_6", t(coef(e_y_35_6)))


e_y_35_7 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_7") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_7)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_35_7", t(coef(e_y_35_7)))


e_y_35_8 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_8") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_8)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_35_8", t(coef(e_y_35_8)))


e_y_70_9 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_70_9") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_9)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_70_9", t(coef(e_y_70_9)))


e_y_70_10 <- aci.prep %>% filter(keep.row == "yes" & 
                                   id == "e_y_70_10" &
                                   Ci < 700) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_10)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_70_10", t(coef(e_y_70_10)))


e_y_70_11 <- aci.prep %>% filter(keep.row == "yes" & 
                                   id == "e_y_70_11" &
                                   Ci < 750) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_11)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_70_11", t(coef(e_y_70_11)))


e_y_70_12 <- aci.prep %>% filter(keep.row == "yes" & 
                                   id == "e_y_70_12" &
                                   Ci < 850) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_12)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_70_12", t(coef(e_y_70_12)))


e_y_105_13 <- aci.prep %>% filter(keep.row == "yes" & 
                                    id == "e_y_105_13" &
                                    Ci < 750) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_13)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_105_13", t(coef(e_y_105_13)))


e_y_105_14 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_105_14") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_14)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_105_14", t(coef(e_y_105_14)))


e_y_105_15 <- aci.prep %>% filter(keep.row == "yes" & 
                                    id == "e_y_105_15" &
                                    Ci < 800) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_15)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_105_15", t(coef(e_y_105_15)))


e_y_105_16 <- aci.prep %>% filter(keep.row == "yes" & 
                                    id == "e_y_105_16" &
                                    Ci < 750) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_16)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_105_16", t(coef(e_y_105_16)))


e_y_140_17 <- aci.prep %>% filter(keep.row == "yes" & 
                                    id == "e_y_140_17" &
                                    Ci < 750) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_17)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_140_17", t(coef(e_y_140_17)))


e_y_140_18 <- aci.prep %>% filter(keep.row == "yes" & 
                                    id == "e_y_140_18" &
                                    Ci < 850) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_18)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_140_18", t(coef(e_y_140_18)))


e_y_140_19 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_140_19") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_19)
aci.coefs[nrow(aci.coefs) + 1,] <- c(id = "e_y_140_19", t(coef(e_y_140_19)))


e_y_140_20 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_140_20") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_20)

e_y_210_21 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_21") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_21)

e_y_210_22 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_22") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_22)

e_y_210_23 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_23") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_23)

e_y_210_24 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_24") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_24)

e_y_280_25 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_25") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_25)

e_y_280_26 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_26") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_26)

e_y_280_27 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_27") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_27)

e_y_280_28 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_28") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_28)

e_y_350_29 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_29") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_29)

e_y_350_30 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_30") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_30)

e_y_350_31 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_31") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_31)

e_y_350_32 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_32") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_32)

e_y_630_33 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_33") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_33)

e_y_630_34 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_34") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_34)

e_y_630_35 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_35") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_35)

e_y_630_36 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_36") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_36)

#######################################
# Elevated CO2 non-inoculated
#######################################
e_n_0_37 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_37") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_37)

e_n_0_38 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_38") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_38)

e_n_0_39 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_39") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_39)

e_n_0_40 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_40") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_40)

e_n_0_40 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_40") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_40)

e_n_35_41 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_41") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_35_41)

e_n_35_42 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_42") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_35_42)

e_n_35_43 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_43") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_35_43)

e_n_35_44 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_44") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_35_44)

e_n_70_45 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_45") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_45)

e_n_70_46 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_46") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_46)

e_n_70_47 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_47") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_47)

e_n_70_48 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_48") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_48)

e_n_105_49 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_49") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_49)

e_n_105_50 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_50") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_50)

e_n_105_51 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_51") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_51)

e_n_105_52 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_52") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_52)

e_n_140_53 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_53") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_53)

e_n_140_54 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_54") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_54)

e_n_140_55 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_55") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_55)

e_n_140_56 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_56") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_56)

e_n_210_57 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_57") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_57)

e_n_210_58 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_58") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_58)

e_n_210_59 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_59") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_59)

e_n_210_60 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_60") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_60)

e_n_280_61 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_61") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_61)

e_n_280_62 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_62") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_62)

e_n_280_63 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_63") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_63)

e_n_280_64 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_64") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_64)

e_n_350_65 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_65") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_350_65)

e_n_350_66 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_66") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_350_66)

e_n_350_67 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_67") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_350_67)

e_n_630_68 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_68") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_68)

e_n_630_69 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_69") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_69)

e_n_630_70 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_70") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_70)

e_n_630_71 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_71") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_71)

e_n_630_72 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_72") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_72)

#######################################
# Ambient CO2 inoculated
#######################################
a_y_0_73 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_73") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
e_n_105_50 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_50") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_50)plot(a_y_0_73)

a_y_0_74 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_74") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_0_74)

a_y_0_75 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_75") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_0_75)

a_y_0_76 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_76") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_0_76)

a_y_35_77 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_77") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_77)

a_y_35_78 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_78") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_78)

a_y_35_79 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_79") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_79)

a_y_35_80 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_80") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_80)

a_y_70_81 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_81") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_81)

a_y_70_82 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_82") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_82)

a_y_70_83 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_83") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_83)

a_y_70_84 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_84") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_84)

a_y_105_85 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_85") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_85)

a_y_105_86 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_86") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_86)

a_y_105_87 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_87") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_87)

a_y_105_88 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_88") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_88)

a_y_140_89 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_89") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_89)

a_y_140_90 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_90") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_90)

a_y_140_91 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_91") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_91)

a_y_140_92 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_92") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_92)

a_y_210_93 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_93") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_ny_210_93)

a_y_210_94 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_94") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_94)

a_y_210_95 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_95") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_95)

a_y_210_96 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_96") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_96)

a_y_280_97 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_97") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_97)

a_y_280_98 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_98") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_98)

a_y_280_99 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_99") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_99)

a_y_280_100 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_100") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_100)

a_y_350_101 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_101") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_101)

a_y_350_102 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_102") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_102)

a_y_350_103 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_103") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_103)

a_y_350_104 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_104") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_104)

a_y_630_105 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_105") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_105)

a_y_630_106 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_106") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_106)

a_y_630_107 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_107") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_107)

a_y_630_108 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_108") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_108)

#######################################
# Ambient CO2 non-inoculated
#######################################
a_n_0_109 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_109") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_0_109)

a_n_0_110 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_110") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_0_110)

a_n_0_111 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_111") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin"),
         fitTPU = TRUE, Tcorrect = FALSE, fitmethod = "bilinear",
         citransition = 500)
plot(a_n_0_111)

a_n_0_112 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_112") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_0_112)

a_n_35_113 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_113") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_113)

a_n_35_114 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_114") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_114)

a_n_35_115 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_115") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_115)

a_n_35_116 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_116") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_116)

a_n_70_117 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_117") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_117)

a_n_70_118 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_118") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_118)

a_n_70_119 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_119") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_119)

a_n_70_120 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_120") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_120)

a_n_105_121 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_121") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_121)

a_n_105_122 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_122") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_122)

a_n_105_123 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_123") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_123)

a_n_105_124 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_124") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_124)

a_n_140_125 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_125") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_125)

a_n_140_126 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_126") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_126)

a_n_140_127 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_127") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_127)

a_n_140_128 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_128") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_128)

a_n_210_129 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_129") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_129)

a_n_210_130 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_130") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_130)

a_n_210_131 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_131") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_131)

a_n_210_132 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_132") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_132)

a_n_280_133 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_133") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_133)

a_n_280_134 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_134") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_134)

a_n_280_135 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_135") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_135)

a_n_280_136 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_136") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_136)

a_n_350_137 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_137") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_137)

a_n_350_138 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_138") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_138)

a_n_350_139 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_139") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_139)

a_n_350_140 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_140") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_140)

a_n_630_141 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_141") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_141)

a_n_630_142 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_142") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_142)

a_n_630_143 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_143") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_143)

a_n_630_144 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_144") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_144)








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