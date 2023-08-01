################################################################
# Load libraries
################################################################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(plantecophys)
library(rpmodel)

################################################################
# Load 
################################################################

################################################################
# Load Lebauer & Treseder (2008) dataset
################################################################
lt.2008 <- read.csv("../data_sheets/lebauer_treseder_table1.csv") %>%
  mutate(biome = factor(biome, 
                        levels = c("wetland", "temp.grassland", "desert",
                                   "trop.grassland", "tundra", "trop.forest",
                                   "temp.forest", "overall")))
head(lt.2008)

################################################################
# Plot LeBauer & Treseder (2008) dataset
################################################################
npp.nlimitation <- ggplot(data = lt.2008, aes(x = eff.size, y = biome)) +
  geom_vline(xintercept = 1, linewidth = 1, color = "red", 
             linetype = "dashed") +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower.ci, xmax = upper.ci), 
                 linewidth = 1, height = 0.5) +
  scale_x_continuous(limits = c(0.5, 2)) +
  scale_y_discrete(labels = c("Wetland", "Temperate grassland", "Desert",
                              "Tropical grassland", "Tundra", "Tropical forest",
                              "Temperate forest", "Overall")) +
  labs(x = "Growth response to nitrogen fertilization",
       y = "Biome") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
npp.nlimitation

png("../working_drafts/figs/NxCO2xI_ESAtalk_lebauer.png",
    height = 4.5, width = 8, res = 600, units = "in")
npp.nlimitation
dev.off()

################################################################
## Make photosynthesis figure (A/Ci for optimality expectation)
################################################################
photo.elv <- data.frame(type = "elv",
                        Photosyn(Ca = seq(50, 2000, 0.25),
                                 Vcmax = 67, Jmax = 150))
photo.amb <- data.frame(type = "amb",
                        Photosyn(Ca = seq(50, 2000, 0.25),
                                 Vcmax = 100, Jmax = 150))

photo.merged <- photo.elv %>% full_join(photo.amb) %>%
  mutate(Anet = pmin(Ac, Aj),
         A_limiting = ifelse(Anet == Ac, "Ac", 
                             ifelse(Anet == Aj, "Aj",
                                    ifelse(Aj == Ac,
                                           "Ac_Aj", NA))))

photo.fig.amb <- ggplot(data = subset(photo.merged, type == "amb"), 
                        aes(x = Ci)) +
  geom_smooth(aes(y = ALEAF), linewidth = 3, color = "black") +
  scale_y_continuous(limits = c(-0.1, 35)) +
  labs(x = expression(bold("Intercellular CO"["2"])),
       y = expression(bold("Net photosynthesis"))) +
  guides(color = "none") +
  theme_bw(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 2))
photo.fig.amb

png("../working_drafts/figs/NxCO2xI_ESAtalk_aci_amb.png",
    width = 6, height = 5, res = 600, units = "in")
photo.fig.amb
dev.off()

photo.fig.both <- ggplot(data = photo.merged, aes(x = Ci)) +
  geom_smooth(aes(y = ALEAF, linetype = type),
              linewidth = 3, color = "black") +
  scale_y_continuous(limits = c(-0.1, 35)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = expression(bold("Intercellular CO"["2"])),
       y = expression(bold("Net photosynthesis"))) +
  guides(color = "none", linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 2))
photo.fig.both

png("../working_drafts/figs/NxCO2xI_ESAtalk_aci_both.png",
    width = 6, height = 5, res = 600, units = "in")
photo.fig.both
dev.off()

################################################################
## P-model
################################################################
co2.sims <- data.frame(co2.acclim = c(420, 1000),
                       anet = c(10, 20),
                       vcmax25 = c(20, 10))
  
  
  rpmodel(tc = 25, vpd = 1.5, co2 = c(420, 1000), 
                    elv = 100, fapar = 0.5, ppfd = 2000) %>%
  data.frame()

vcmax.optimal <- ggplot(data = co2.sims, aes(x = as.factor(co2.acclim), 
                                             y = vcmax25)) +
  geom_point(size = 4) +
  scale_y_continuous(limits = c(7.5, 22.5)) +
  labs(x = expression(bold("Acclimated CO"[2]*" concentration (ppm)")),
       y = "Photosynthetic capacity") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"))

anet.optimal <- ggplot(data = co2.sims, aes(x = as.factor(co2.acclim), 
                                           y = anet)) +
  geom_point(size = 4) +
  scale_y_continuous(limits = c(7.5, 22.5)) +
  labs(x = expression(bold("Acclimated CO"[2]*" concentration (ppm)")),
       y = "Net photosynthesis") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"))

png("../working_drafts/figs/NxCO2xI_ESAtalk_optimality_preds.png",
    height = 4.5, width = 10, res = 600, units = "in")
ggarrange(anet.optimal, vcmax.optimal, ncol = 2, align = "hv")
dev.off()
