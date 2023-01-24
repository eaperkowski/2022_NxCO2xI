###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
co2.response.wk6 <- read.csv("../data_sheets/NxCO2_co2_resp_wk6.csv")
rd.wk6 <- read.csv("../data_sheets/NxCO2_rd_wk7.csv")

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
aci.prep <- co2.response.wk6 %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair,	Tleaf) %>%
  arrange(id) %>%
  left_join(rd.wk6, by = "id") %>%
  dplyr::select(-week, -tLeaf) %>%
  group_by(id) %>%
  mutate(keep.row = "yes") %>%
  data.frame()

aci.temps <- aci.prep %>%
  group_by(id) %>%
  summarize(Tleaf = mean(Tleaf, na.rm = TRUE))

#####################################################################
# A/Ci curves
#####################################################################

#######################################
# Elevated CO2 inoculated
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


e_y_0_2 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_0_2")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_2)
aci.coefs[2,] <- c(id = "e_y_0_2", t(coef(e_y_0_2)))


e_y_0_3 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_0_3") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_3)
aci.coefs[3,] <- c(id = "e_y_0_3", t(coef(e_y_0_3)))


e_y_0_4 <- aci.prep %>% filter(keep.row == "yes" & 
                                 id == "e_y_0_4") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_0_4)
aci.coefs[4,] <- c(id = "e_y_0_4", t(coef(e_y_0_4)))


e_y_35_5 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_5") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_5)
aci.coefs[5,] <- c(id = "e_y_35_5", t(coef(e_y_35_5)))


e_y_35_6 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_6") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_6)
aci.coefs[6,] <- c(id = "e_y_35_6", t(coef(e_y_35_6)))


e_y_35_7 <- aci.prep %>% filter(keep.row == "yes" & 
                                  id == "e_y_35_7" & CO2_r < 870)  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_7)
aci.coefs[7,] <- c(id = "e_y_35_7", t(coef(e_y_35_7)))


e_y_35_8 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_35_8")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_35_8)
aci.coefs[8,] <- c(id = "e_y_35_8", t(coef(e_y_35_8)))


e_y_70_9 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_70_9") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_9)
aci.coefs[9,] <- c(id = "e_y_70_9", t(coef(e_y_70_9)))


e_y_70_10 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_70_10") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_10)
aci.coefs[10,] <- c(id = "e_y_70_10", t(coef(e_y_70_10)))


e_y_70_11 <- aci.prep %>% filter(keep.row == "yes" & 
                                   id == "e_y_70_11") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_11)
aci.coefs[11,] <- c(id = "e_y_70_11", t(coef(e_y_70_11)))


e_y_70_12 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_70_12") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_70_12)
aci.coefs[12,] <- c(id = "e_y_70_12", t(coef(e_y_70_12)))


e_y_105_13 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_105_13")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_13)
aci.coefs[13,] <- c(id = "e_y_105_13", t(coef(e_y_105_13)))


e_y_105_14 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_105_14") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_14)
aci.coefs[14,] <- c(id = "e_y_105_14", t(coef(e_y_105_14)))


e_y_105_15 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_105_15") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_15)
aci.coefs[15,] <- c(id = "e_y_105_15", t(coef(e_y_105_15)))


e_y_105_16 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_105_16") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_105_16)
aci.coefs[16,] <- c(id = "e_y_105_16", t(coef(e_y_105_16)))


e_y_140_17 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_140_17") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_17)
aci.coefs[17,] <- c(id = "e_y_140_17", t(coef(e_y_140_17)))


e_y_140_18 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_140_18") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_18)
aci.coefs[18,] <- c(id = "e_y_140_18", t(coef(e_y_140_18)))


e_y_140_19 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_140_19") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_19)
aci.coefs[19,] <- c(id = "e_y_140_19", t(coef(e_y_140_19)))


e_y_140_20 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_140_20") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_140_20)
aci.coefs[20,] <- c(id = "e_y_140_20", t(coef(e_y_140_20)))


e_y_210_21 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_21") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_21)
aci.coefs[21,] <- c(id = "e_y_210_21", t(coef(e_y_210_21)))


e_y_210_22 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_22" & Ci < 800) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_22)
aci.coefs[22,] <- c(id = "e_y_210_22", t(coef(e_y_210_22)))


e_y_210_23 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_23") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_23)
aci.coefs[23,] <- c(id = "e_y_210_23", t(coef(e_y_210_23)))


e_y_210_24 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_210_24") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_210_24)
aci.coefs[24,] <- c(id = "e_y_210_24", t(coef(e_y_210_24)))


e_y_280_25 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_25") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_25)
aci.coefs[25,] <- c(id = "e_y_280_25", t(coef(e_y_280_25)))


e_y_280_26 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_26") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_26)
aci.coefs[26,] <- c(id = "e_y_280_26", t(coef(e_y_280_26)))


e_y_280_27 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_27" & Ci < 800) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_280_27)
aci.coefs[27,] <- c(id = "e_y_280_27", t(coef(e_y_280_27)))


e_y_280_28 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_280_28") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_y_280_28)
aci.coefs[28,] <- c(id = "e_y_280_28", t(coef(e_y_280_28)))


e_y_350_29 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_29") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_29)
aci.coefs[29,] <- c(id = "e_y_350_29", t(coef(e_y_350_29)))


e_y_350_30 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_30") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_30)
aci.coefs[30,] <- c(id = "e_y_350_30", t(coef(e_y_350_30)))


e_y_350_31 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_31")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_31)
aci.coefs[31,] <- c(id = "e_y_350_31", t(coef(e_y_350_31)))


e_y_350_32 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_350_32")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_350_32)
aci.coefs[32,] <- c(id = "e_y_350_32", t(coef(e_y_350_32)))

e_y_630_33 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_33") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_33)
aci.coefs[33,] <- c(id = "e_y_630_33", t(coef(e_y_630_33)))

e_y_630_34 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_34") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_34)
aci.coefs[34,] <- c(id = "e_y_630_34", t(coef(e_y_630_34)))


e_y_630_35 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_35") %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_35)
aci.coefs[35,] <- c(id = "e_y_630_35", t(coef(e_y_630_35)))

e_y_630_36 <- aci.prep %>% filter(keep.row == "yes" & id == "e_y_630_36") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_y_630_36)
aci.coefs[36,] <- c(id = "e_y_630_36", t(coef(e_y_630_36)))

#######################################
# Elevated CO2 non-inoculated
#######################################
e_n_0_37 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_37" & A > 4) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin",
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_37)
aci.coefs[37,] <- c(id = "e_n_0_37", t(coef(e_n_0_37)))


e_n_0_38 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_38") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin",
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_38)
aci.coefs[38,] <- c(id = "e_n_0_38", t(coef(e_n_0_38)))


e_n_0_39 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_39") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin",
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_0_39)
aci.coefs[39,] <- c(id = "e_n_0_39", t(coef(e_n_0_39)))

e_n_0_40 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_0_40" & A > 2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_0_40)
aci.coefs[40,] <- c(id = "e_n_0_40", t(coef(e_n_0_40)))

e_n_35_41 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_41")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_35_41)
aci.coefs[41,] <- c(id = "e_n_35_41", t(coef(e_n_35_41)))

e_n_35_42 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_42") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin",
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_35_42)
aci.coefs[42,] <- c(id = "e_n_35_42", t(coef(e_n_35_42)))


e_n_35_43 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_43") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_35_43)
aci.coefs[43,] <- c(id = "e_n_35_43", t(coef(e_n_35_43)))

e_n_35_44 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_35_44") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_35_44)
aci.coefs[44,] <- c(id = "e_n_35_44", t(coef(e_n_35_44)))

e_n_70_45 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_45")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_45)
aci.coefs[45,] <- c(id = "e_n_70_45", t(coef(e_n_70_45)))

e_n_70_46 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_46") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_70_46)
aci.coefs[46,] <- c(id = "e_n_70_46", t(coef(e_n_70_46)))

e_n_70_47 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_47") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_47)
aci.coefs[47,] <- c(id = "e_n_70_47", t(coef(e_n_70_47)))

e_n_70_48 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_70_48") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_70_48)
aci.coefs[48,] <- c(id = "e_n_70_48", t(coef(e_n_70_48)))


e_n_105_49 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_49") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_105_49)
aci.coefs[49,] <- c(id = "e_n_105_49", t(coef(e_n_105_49)))


e_n_105_50 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_50") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_50)
aci.coefs[50,] <- c(id = "e_n_105_50", NA, NA, NA, NA)


e_n_105_51 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_61") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_51)
aci.coefs[51,] <- c(id = "e_n_105_51", t(coef(e_n_105_51)))


e_n_105_52 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_105_52") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_105_52)
aci.coefs[52,] <- c(id = "e_n_105_52", t(coef(e_n_105_52)))


e_n_140_53 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_53") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_53)
aci.coefs[53,] <- c(id = "e_n_140_53", t(coef(e_n_140_53)))


e_n_140_54 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_54" & Ci < 700) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_54)
aci.coefs[54,] <- c(id = "e_n_140_54", t(coef(e_n_140_54)))

e_n_140_55 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_55") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_55)
aci.coefs[55,] <- c(id = "e_n_140_55", t(coef(e_n_140_55)))

e_n_140_56 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_140_56") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_140_56)
aci.coefs[56,] <- c(id = "e_n_140_56", t(coef(e_n_140_56)))


e_n_210_57 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_57") %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_57)
aci.coefs[57,] <- c(id = "e_n_210_57", t(coef(e_n_210_57)))

e_n_210_58 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_58") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_210_58)
aci.coefs[58,] <- c(id = "e_n_210_58", t(coef(e_n_210_58)))


e_n_210_59 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_59") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_59)
aci.coefs[59,] <- c(id = "e_n_210_59", t(coef(e_n_210_59)))


e_n_210_60 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_210_60") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_210_60)
aci.coefs[60,] <- c(id = "e_n_210_60", t(coef(e_n_210_60)))


e_n_280_61 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_61") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_61)
aci.coefs[61,] <- c(id = "e_n_280_61", t(coef(e_n_280_61)))


e_n_280_62 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_62" & A > 2.5) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_62)
aci.coefs[62,] <- c(id = "e_n_280_62", t(coef(e_n_280_62)))


e_n_280_63 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_63") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_280_63)
aci.coefs[63,] <- c(id = "e_n_280_63", t(coef(e_n_280_63)))


e_n_280_64 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_280_64") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_280_64)
aci.coefs[64,] <- c(id = "e_n_280_64", t(coef(e_n_280_64)))

e_n_350_65 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_65")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_350_65)
aci.coefs[65,] <- c(id = "e_n_350_65", t(coef(e_n_350_65)))


e_n_350_66 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_66")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_350_66)
aci.coefs[66,] <- c(id = "e_n_350_66", t(coef(e_n_350_66)))

e_n_350_67 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_67") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(e_n_350_67)
aci.coefs[67,] <- c(id = "e_n_350_67", t(coef(e_n_350_67)))


e_n_350_68 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_350_68") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_350_68)
aci.coefs[68,] <- c(id = "e_n_350_68", t(coef(e_n_350_68)))


e_n_630_69 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_69") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_69)
aci.coefs[69,] <- c(id = "e_n_630_69", t(coef(e_n_630_69)))


e_n_630_70 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_70")  %>%
  mutate(A = A/2) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_70)
aci.coefs[70,] <- c(id = "e_n_630_70", t(coef(e_n_630_70)))


e_n_630_71 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_71") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_71)
aci.coefs[71,] <- c(id = "e_n_630_71", t(coef(e_n_630_71)))


e_n_630_72 <- aci.prep %>% filter(keep.row == "yes" & id == "e_n_630_72") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(e_n_630_72)
aci.coefs[72,] <- c(id = "e_n_630_72", t(coef(e_n_630_72)))

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
plot(a_y_0_73)
aci.coefs[73,] <- c(id = "a_y_0_73", t(coef(a_y_0_73)))


a_y_0_74 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_74") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_0_74)
aci.coefs[74,] <- c(id = "a_y_0_74", t(coef(a_y_0_74)))

a_y_0_75 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_75") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_0_75)
aci.coefs[75,] <- c(id = "a_y_0_75", t(coef(a_y_0_75)))

a_y_0_76 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_0_76") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_0_76)
aci.coefs[76,] <- c(id = "a_y_0_76", t(coef(a_y_0_76)))


a_y_35_77 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_77") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_77)
aci.coefs[77,] <- c(id = "a_y_35_77", t(coef(a_y_35_77)))


a_y_35_78 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_78") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_78)
aci.coefs[78,] <- c(id = "a_y_35_78", t(coef(a_y_35_78)))


a_y_35_79 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_79") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_79)
aci.coefs[79,] <- c(id = "a_y_35_79", t(coef(a_y_35_78)))


a_y_35_80 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_35_80") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_35_80)
aci.coefs[80,] <- c(id = "a_y_35_80", t(coef(a_y_35_80)))


a_y_70_81 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_81") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_81)
aci.coefs[81,] <- c(id = "a_y_70_81", t(coef(a_y_70_81)))


a_y_70_82 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_82") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_82)
aci.coefs[82,] <- c(id = "a_y_70_82", t(coef(a_y_70_82)))

a_y_70_83 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_83") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_83)
aci.coefs[83,] <- c(id = "a_y_70_83", t(coef(a_y_70_83)))

a_y_70_84 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_70_84") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_70_84)
aci.coefs[84,] <- c(id = "a_y_70_84", t(coef(a_y_70_84)))

a_y_105_85 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_85" & A < 40) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_85)
aci.coefs[85,] <- c(id = "a_y_105_85", t(coef(a_y_105_85)))


a_y_105_86 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_86") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_86)
aci.coefs[86,] <- c(id = "a_y_105_86", t(coef(a_y_105_86)))


a_y_105_87 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_87") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_87)
aci.coefs[87,] <- c(id = "a_y_105_87", t(coef(a_y_105_87)))


a_y_105_88 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_105_88") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_105_88)
aci.coefs[88,] <- c(id = "a_y_105_88", t(coef(a_y_105_88)))

a_y_140_89 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_89") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_89)
aci.coefs[89,] <- c(id = "a_y_140_89", t(coef(a_y_140_89)))


a_y_140_90 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_90") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_90)
aci.coefs[90,] <- c(id = "a_y_140_90", t(coef(a_y_140_90)))


a_y_140_91 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_91") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_91)
aci.coefs[91,] <- c(id = "a_y_140_91", t(coef(a_y_140_91)))

a_y_140_92 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_140_92") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_140_92)
aci.coefs[92,] <- c(id = "a_y_140_92", t(coef(a_y_140_92)))

a_y_210_93 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_93") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_93)
aci.coefs[93,] <- c(id = "a_y_210_93", t(coef(a_y_210_93)))


a_y_210_94 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_94") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_94)
aci.coefs[94,] <- c(id = "a_y_210_94", t(coef(a_y_210_94)))


a_y_210_95 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_95") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_95)
aci.coefs[95,] <- c(id = "a_y_210_95", t(coef(a_y_210_95)))


a_y_210_96 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_210_96") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_210_96)
aci.coefs[96,] <- c(id = "a_y_210_96", t(coef(a_y_210_96)))


a_y_280_97 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_97") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_97)
aci.coefs[97,] <- c(id = "a_y_280_97", t(coef(a_y_280_97)))


a_y_280_98 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_98") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_98)
aci.coefs[98,] <- c(id = "a_y_280_98", t(coef(a_y_280_98)))


a_y_280_99 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_99") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_99)
aci.coefs[99,] <- c(id = "a_y_280_99", t(coef(a_y_280_99)))


a_y_280_100 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_280_100") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_280_100)
aci.coefs[100,] <- c(id = "a_y_280_100", t(coef(a_y_280_100)))


a_y_350_101 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_101") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_101)
aci.coefs[101,] <- c(id = "a_y_350_101", t(coef(a_y_350_101)))


a_y_350_102 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_102") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_102)
aci.coefs[102,] <- c(id = "a_y_350_102", t(coef(a_y_350_102)))


a_y_350_103 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_103") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_103)
aci.coefs[103,] <- c(id = "a_y_350_103", t(coef(a_y_350_103)))


a_y_350_104 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_350_104") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_350_104)
aci.coefs[104,] <- c(id = "a_y_350_104", t(coef(a_y_350_104)))


a_y_630_105 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_105") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_105)
aci.coefs[105,] <- c(id = "a_y_630_105", t(coef(a_y_630_105)))


a_y_630_106 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_106") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_106)
aci.coefs[106,] <- c(id = "a_y_630_106", t(coef(a_y_630_106)))


a_y_630_107 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_107") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_107)
aci.coefs[107,] <- c(id = "a_y_630_107", t(coef(a_y_630_107)))


a_y_630_108 <- aci.prep %>% filter(keep.row == "yes" & id == "a_y_630_108") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_y_630_108)
aci.coefs[108,] <- c(id = "a_y_630_108", t(coef(a_y_630_108)))


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
aci.coefs[109,] <- c(id = "a_n_0_109", t(coef(a_n_0_109)))


a_n_0_110 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_110") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_0_110)
aci.coefs[110,] <- c(id = "a_n_0_110", t(coef(a_n_0_110)))


a_n_0_111 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_111") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin",
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(a_n_0_111)
aci.coefs[111,] <- c(id = "a_n_0_111", t(coef(a_n_0_111)))


a_n_0_112 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_0_112") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_0_112)
aci.coefs[112,] <- c(id = "a_n_0_112", t(coef(a_n_0_112)))


a_n_35_113 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_113") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_113)
aci.coefs[113,] <- c(id = "a_n_35_113", NA, NA, NA, NA)


a_n_35_114 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_114") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_114)
aci.coefs[114,] <- c(id = "a_n_35_114", t(coef(a_n_35_114)))


a_n_35_115 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_115") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_115)
aci.coefs[115,] <- c(id = "a_n_35_115", NA, NA, NA, NA)


a_n_35_116 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_35_116") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_35_116)
aci.coefs[116,] <- c(id = "a_n_35_116", t(coef(a_n_35_116)))


a_n_70_117 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_117") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_117)
aci.coefs[117,] <- c(id = "a_n_70_117", t(coef(a_n_70_117)))


a_n_70_118 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_118") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_118)
aci.coefs[118,] <- c(id = "a_n_70_118", t(coef(a_n_70_118)))


a_n_70_119 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_119") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_119)
aci.coefs[119,] <- c(id = "a_n_70_119", t(coef(a_n_70_119)))


a_n_70_120 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_70_120") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_70_120)
aci.coefs[120,] <- c(id = "a_n_70_120", t(coef(a_n_70_120)))


a_n_105_121 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_121") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_121)
aci.coefs[121,] <- c(id = "a_n_105_121", t(coef(a_n_105_121)))


a_n_105_122 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_122") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_122)
aci.coefs[122,] <- c(id = "a_n_105_122", t(coef(a_n_105_122)))


a_n_105_123 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_123") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_123)
aci.coefs[123,] <- c(id = "a_n_105_123", t(coef(a_n_105_123)))


a_n_105_124 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_105_124") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_105_124)
aci.coefs[124,] <- c(id = "a_n_105_124", t(coef(a_n_105_124)))


a_n_140_125 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_125") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_125)
aci.coefs[125,] <- c(id = "a_n_140_125", t(coef(a_n_140_125)))


a_n_140_126 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_126" & A < 40) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_126)
aci.coefs[126,] <- c(id = "a_n_140_126", t(coef(a_n_140_126)))


a_n_140_127 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_127") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_127)
aci.coefs[127,] <- c(id = "a_n_140_127", t(coef(a_n_140_127)))


a_n_140_128 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_140_128") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_140_128)
aci.coefs[128,] <- c(id = "a_n_140_128", t(coef(a_n_140_128)))


a_n_210_129 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_129") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_129)
aci.coefs[129,] <- c(id = "a_n_210_129", t(coef(a_n_210_129)))


a_n_210_130 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_130") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_130)
aci.coefs[130,] <- c(id = "a_n_210_130", t(coef(a_n_210_130)))


a_n_210_131 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_131") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_131)
aci.coefs[131,] <- c(id = "a_n_210_131", t(coef(a_n_210_131)))


a_n_210_132 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_210_132") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_210_132)
aci.coefs[132,] <- c(id = "a_n_210_132", t(coef(a_n_210_132)))


a_n_280_133 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_133") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_133)
aci.coefs[133,] <- c(id = "a_n_280_133", t(coef(a_n_280_133)))


a_n_280_134 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_134") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_134)
aci.coefs[134,] <- c(id = "a_n_280_134", t(coef(a_n_280_134)))


a_n_280_135 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_135") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_135)
aci.coefs[135,] <- c(id = "a_n_280_135", t(coef(a_n_280_135)))


a_n_280_136 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_280_136") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_280_136)
aci.coefs[136,] <- c(id = "a_n_280_136", t(coef(a_n_280_136)))


a_n_350_137 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_137") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_137)
aci.coefs[137,] <- c(id = "a_n_350_137", t(coef(a_n_350_137)))


a_n_350_138 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_138") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_138)
aci.coefs[138,] <- c(id = "a_n_350_138", t(coef(a_n_350_138)))


a_n_350_139 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_139") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_139)
aci.coefs[139,] <- c(id = "a_n_350_139", t(coef(a_n_350_139)))


a_n_350_140 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_350_140") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_350_140)
aci.coefs[140,] <- c(id = "a_n_350_140", t(coef(a_n_350_140)))


a_n_630_141 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_141") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_141)
aci.coefs[141,] <- c(id = "a_n_630_141", t(coef(a_n_630_141)))


a_n_630_142 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_142") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_142)
aci.coefs[142,] <- c(id = "a_n_630_142", t(coef(a_n_630_142)))


a_n_630_143 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_143") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_143)
aci.coefs[143,] <- c(id = "a_n_630_143", t(coef(a_n_630_143)))


a_n_630_144 <- aci.prep %>% filter(keep.row == "yes" & id == "a_n_630_144") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(a_n_630_144)
aci.coefs[144,] <- c(id = "a_n_630_144", t(coef(a_n_630_144)))


#####################################################################
# A/Ci curves
#####################################################################
aci.coefs$Vcmax <- as.numeric(aci.coefs$Vcmax)
aci.coefs$Jmax <- as.numeric(aci.coefs$Jmax)
aci.coefs$Rd <- as.numeric(aci.coefs$Rd)
aci.coefs$TPU <- as.numeric(aci.coefs$TPU)

aci.coefs[, c(2:5)] <- round(aci.coefs[, c(2:5)], digits = 3)

aci.fits <- aci.coefs %>% left_join(aci.temps) %>%
  mutate(Jmax = ifelse(Jmax > 300, NA, Jmax),
         vcmax25 = temp_standardize(Vcmax, "Vcmax", standard.to = 25,
                                    tLeaf = Tleaf, tGrow = 22.5),
         jmax25 = temp_standardize(Jmax, "Jmax", standard.to = 25,
                                   tLeaf = Tleaf, tGrow = 22.5),
         jmax25.vcmax25 = jmax25 / vcmax25) %>%
  select(id, tleaf = Tleaf, vcmax25, jmax25, jmax25.vcmax25, rd25 = Rd, tpu = TPU) %>%
  mutate_if(is.numeric, round, 3)

anet <- aci.prep %>%
  filter(CO2_r > 419.5 & CO2_r < 420.5) %>%
  group_by(id) %>%
  summarize(anet = mean(A),
            gsw = mean(gsw),
            ci.ca = mean(Ci) / mean(Ca)) %>%
  mutate(id = ifelse(id == "e_n_105_61", "e_n_105_51", id))


photo.data <- anet %>% full_join(aci.fits) %>% 
  mutate_if(is.numeric, round, 3) %>%
  mutate(jmax25 = ifelse(jmax25 > 300, NA, jmax25),
         jmax25) %>%
  select(id, tleaf, anet,  gsw, ci.ca, vcmax25, jmax25, 
         jmax25.vcmax25, rd25, tpu) %>%
  filter (id != "e_y_35_7_b") %>%
  setNames(c(names(.)[1], paste0(names(.)[-1],"_wk6")))

write.csv(photo.data, "../data_sheets/NxCO2_photo_data_wk6.csv", row.names = FALSE)

full.df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  full_join(photo.data)

write.csv(full.df, "../data_sheets/NxCO2xI_compiled_datasheet.csv", row.names = FALSE)




