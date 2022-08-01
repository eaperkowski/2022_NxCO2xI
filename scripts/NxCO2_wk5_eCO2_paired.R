## Load libraries
library(plantecophys)
library(readLicorData)
library(ggplot2)
library(dplyr)
library(stringr)

## Load and clean raw licor file
licor.file <- licorData("/Users/eaperkowski/git/2022_NxCO2_growthChamber/licor_raw/week5/2022-07-29_paired_W5D3")
licor.file$id[193:288] <- "e_n_630_71_DAT"

## Clean licor file cols and reorganize cols. Also separate id and curve type
## into separate cols
licor.file <- as.data.frame(licor.file) %>%
  select(obs, date, curve.id = id, A, Adyn, Asty, Ca, Ci, gsw, VPDleaf, Qin, Tleaf) %>%
  mutate(curve.type = str_extract(curve.id, "[^_]+$"),
         id = gsub("_[^_]+$", "", curve.id))

## Visualize curves
e_y_630 <- ggplot(data = subset(licor.file, id == "e_y_630_35"), 
       aes(x = Ci, y = A, fill = curve.type)) +
  geom_point(shape = 21) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_x_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) +
  theme_bw(base_size = 18) +
  labs(fill = "Curve type",
       title = "E_Y_630_35")

e_n_630 <- ggplot(data = subset(licor.file, id == "e_n_630_71"), 
       aes(x = Ci, y = A, fill = curve.type)) +
  geom_point(shape = 21) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_x_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) +
  theme_bw(base_size = 18) +
  labs(fill = "Curve type",
       title = "E_N_630_71")

e_y_0 <- ggplot(data = subset(licor.file, id == "e_y_0_2"), 
       aes(x = Ci, y = A, fill = curve.type)) +
  geom_point(shape = 21) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_x_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) +
  theme_bw(base_size = 18) +
  labs(fill = "Curve type",
       title = "E_Y_0_2")

e_y_210 <- ggplot(data = subset(licor.file, id == "e_y_210_24"), 
       aes(x = Ci, y = A, fill = curve.type)) +
  geom_point(shape = 21) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_x_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) +
  theme_bw(base_size = 18) +
  labs(fill = "Curve type",
       title = "E_Y_210_24")

e_n_210 <- ggplot(data = subset(licor.file, id == "e_n_210_59"), 
                  aes(x = Ci, y = A, fill = curve.type)) +
  geom_point(shape = 21) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_x_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) +
  theme_bw(base_size = 18) +
  labs(fill = "Curve type",
       title = "E_N_210_59")

## Write figs to central file
png("../docs/NxCO2_wk5_eCO2_paired.png",
    width = 12, height = 12, units = 'in', res = 600)
ggpubr::ggarrange(e_n_630, e_y_630, e_n_210, e_y_210, e_y_0, 
                  legend = "right", common.legend = TRUE, ncol = 2, nrow = 3)
dev.off()

## Let's do some curve fitting! Note: curves are temp standardized using built-in
## fxns (will be manually standardized for actual model code)

# E_Y_0
e_y_0_DATfit <- licor.file %>% filter(curve.id == "e_y_0_2_DAT" & Ca < 700) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_y_0_DATfit)


summary(e_y_0_DATfit)[[]]

e_y_0_dat_results <- data.frame(id = "e_y_0_2", 
                                type = "DAT",
                                vcmax25 = e_y_0_DATfit$pars[1],
                                vcmax25_se = e_y_0_DATfit$pars[4],
                                jmax25 = e_y_0_DATfit$pars[2],
                                jmax25_se = e_y_0_DATfit$pars[5],
                                row.names = NULL)

e_y_0_SSfit <- licor.file %>% filter(curve.id == "e_y_0_2_steady" & Ca < 800) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_y_0_SSfit)

e_y_0_SS_results <- data.frame(id = "e_y_0_2", 
                               type = "steady",
                               vcmax25 = e_y_0_SSfit$pars[1],
                               vcmax25_se = e_y_0_SSfit$pars[4],
                               jmax25 = e_y_0_SSfit$pars[2],
                               jmax25_se = e_y_0_SSfit$pars[5],
                               row.names = NULL)

# E_Y_210
e_y_210_DATfit <- licor.file %>% filter(curve.id == "e_y_210_24_DAT" & Ci < 600) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_y_210_DATfit)

e_y_210_dat_results <- data.frame(id = "e_y_210_24", 
                                  type = "DAT",
                                  vcmax25 = e_y_210_DATfit$pars[1],
                                  vcmax25_se = e_y_210_DATfit$pars[4],
                                  jmax25 = e_y_210_DATfit$pars[2],
                                  jmax25_se = e_y_210_DATfit$pars[5],
                                  row.names = NULL)

e_y_210_SSfit <- licor.file %>% filter(curve.id == "e_y_210_24_steady" & Ci < 800) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
         fitTPU = TRUE)
plot(e_y_0_SSfit)

e_y_210_SS_results <- data.frame(id = "e_y_210_24",
                                 type = "steady",
                                 vcmax25 = e_y_210_SSfit$pars[1],
                                 vcmax25_se = e_y_210_SSfit$pars[4],
                                 jmax25 = e_y_210_SSfit$pars[2],
                                 jmax25_se = e_y_210_SSfit$pars[5],
                                 row.names = NULL)

# E_N_210
e_n_210_DATfit <- licor.file %>% filter(curve.id == "e_n_210_59_DAT" & Ci < 600) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_n_210_DATfit)

e_n_210_dat_results <- data.frame(id = "e_n_210_59", 
                                  type = "DAT",
                                  vcmax25 = e_n_210_DATfit$pars[1],
                                  vcmax25_se = e_n_210_DATfit$pars[4],
                                  jmax25 = e_n_210_DATfit$pars[2],
                                  jmax25_se = e_n_210_DATfit$pars[5],
                                  row.names = NULL)

e_n_210_SSfit <- licor.file %>% filter(curve.id == "e_n_210_59_steady" & Ci < 600) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_n_210_SSfit)

e_n_210_SS_results <- data.frame(id = "e_n_210_59",
                                 type = "steady",
                                 vcmax25 = e_n_210_SSfit$pars[1],
                                 vcmax25_se = e_n_210_SSfit$pars[4],
                                 jmax25 = e_n_210_SSfit$pars[2],
                                 jmax25_se = e_n_210_SSfit$pars[5],
                                 row.names = NULL)

# E_N_630
e_n_630_DATfit <- licor.file %>% filter(curve.id == "e_n_630_71_DAT" & Ci < 600) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_n_630_DATfit)

e_n_630_dat_results <- data.frame(id = "e_n_630_71",
                                  type = "DAT",
                                  vcmax25 = e_n_630_DATfit$pars[1],
                                  vcmax25_se = e_n_630_DATfit$pars[4],
                                  jmax25 = e_n_630_DATfit$pars[2],
                                  jmax25_se = e_n_630_DATfit$pars[5],
                                  row.names = NULL)

e_n_630_SSfit <- licor.file %>% filter(curve.id == "e_n_630_71_steady" & Ci < 800) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_n_630_SSfit)

e_n_630_SS_results <- data.frame(id = "e_n_630_71",
                                 type = "steady",
                                 vcmax25 = e_n_630_SSfit$pars[1],
                                 vcmax25_se = e_n_630_SSfit$pars[4],
                                 jmax25 = e_n_630_SSfit$pars[2],
                                 jmax25_se = e_n_630_SSfit$pars[5],
                                 row.names = NULL)

# E_Y_630
e_y_630_DATfit <- licor.file %>% filter(curve.id == "e_y_630_35_DAT" & Ci < 700) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_y_630_DATfit)

e_y_630_dat_results <- data.frame(id = "e_y_630_35",
                                  type = "DAT",
                                  vcmax25 = e_y_630_DATfit$pars[1],
                                  vcmax25_se = e_y_630_DATfit$pars[4],
                                  jmax25 = e_y_630_DATfit$pars[2],
                                  jmax25_se = e_y_630_DATfit$pars[5],
                                  row.names = NULL)

e_y_630_SSfit <- licor.file %>% filter(curve.id == "e_y_630_35_steady" & Ci < 900) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"))
plot(e_y_630_SSfit)

e_y_630_SS_results <- data.frame(id = "e_y_630_35",
                                 type = "steady",
                                 vcmax25 = e_y_630_SSfit$pars[1],
                                 vcmax25_se = e_y_630_SSfit$pars[4],
                                 jmax25 = e_y_630_SSfit$pars[2],
                                 jmax25_se = e_y_630_SSfit$pars[5],
                                 row.names = NULL)

## Merge DAT and ss dfs
fit_results <- e_y_0_dat_results %>%
  full_join(e_y_0_SS_results) %>%
  full_join(e_y_210_dat_results) %>%
  full_join(e_y_210_SS_results) %>%
  full_join(e_n_210_dat_results) %>%
  full_join(e_n_210_SS_results) %>%
  full_join(e_y_630_dat_results) %>%
  full_join(e_y_630_SS_results) %>%
  full_join(e_n_630_dat_results) %>%
  full_join(e_n_630_SS_results) %>%
  mutate(id = factor(id, levels = c("e_y_0_2", "e_n_210_59", "e_y_210_24",
                                    "e_n_630_71", "e_y_630_35")))


vcmax_comps <- ggplot(data = fit_results, aes(x = id, fill = type)) +
  geom_bar(aes(y = vcmax25), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = vcmax25 - vcmax25_se, ymax = vcmax25 + vcmax25_se, 
                    group = type), position = position_dodge(width = 0.9), 
                width = 0.4, size = 1) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("e_y_0", "e_n_210", "e_y_210",
                              "e_n_630", "e_y_630")) +
  labs(x = "Pot ID", y = expression("V"[cmax25]~ "(μmol m"^-2~"s"^-1~")"),
       fill = "Curve type") +
  theme_bw(base_size = 18)

jmax_comps <- ggplot(data = fit_results, aes(x = id, fill = type)) +
  geom_bar(aes(y = jmax25), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = jmax25 - jmax25_se, ymax = jmax25 + jmax25_se, 
                    group = type), position = position_dodge(width = 0.9), 
                width = 0.4, size = 1) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("e_y_0", "e_n_210", "e_y_210",
                              "e_n_630", "e_y_630")) +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
  labs(x = "Pot ID", y = expression("J"[max25]~ "(μmol m"^-2~"s"^-1~")"),
       fill = "Curve type") +
  theme_bw(base_size = 18)

png("../docs/NxCO2_wk5_eCO2_paired_rateEstimates.png",
    width = 16, height = 8, units = 'in', res = 600)

ggpubr::ggarrange(vcmax_comps, jmax_comps, ncol = 2, common.legend = TRUE,
                  legend = "right")
dev.off()
