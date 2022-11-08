## Load libraries
library(ggplot2)
library(car)
library(emmeans)
library(ggpubr)

## Load file
df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  mutate(n.trt = as.numeric(n.trt),
         co2 = ifelse(co2 == "a", "amb", "elv"),
         inoc = ifelse(inoc == "n", "no.inoc", "inoc"))
  
######################################################################
## Evidence against rolling rest of chlorophyll leaves
######################################################################
narea.comp <- ggplot(data = df, aes(x = narea, y = narea.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc)) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_cor(label.y = 3) +
  stat_regline_equation(label.y = 2.8) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("inoculated",
                                "not inoculated")) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("ambient", "elevated")) +
  scale_x_continuous(limits = c(0,3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(0,3), breaks = seq(0, 3, 1)) +
  labs(x = expression("N"["area"]*" (gN m"^"-2"*")"),
       y = expression("N"["area_chl"]*" (gN m"^"-2"*")"),
       fill = expression("CO"[2]*" treatment"),
       shape = "Inoculation") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw(base_size = 20)

marea.comp <- ggplot(data = subset(df, marea.chl < 100), 
       aes(x = marea, y = marea.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc)) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_cor(label.y = 100) +
  stat_regline_equation(label.y = 95) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("inoculated",
                                "not inoculated")) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("ambient", "elevated")) +
  scale_x_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  scale_y_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  labs(x = expression("M"["area"]*" (g m"^"-2"*")"),
       y = expression("M"["area_chl"]*" (g m"^"-2"*")"),
       fill = expression("CO"[2]*" treatment"),
       shape = "Inoculation") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw(base_size = 18)

nmass.comp <- ggplot(data = df, aes(x = nmass.focal, y = nmass.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc)) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_cor(label.y = 0.075) +
  stat_regline_equation(label.y = 0.070) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("inoculated",
                                "not inoculated")) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("ambient", "elevated")) +
  scale_x_continuous(limits = c(0, 0.075), breaks = seq(0, 0.075, 0.025)) +
  scale_y_continuous(limits = c(0, 0.075), breaks = seq(0, 0.075, 0.025)) +
  labs(x = expression("N"["mass"]*" (gN m"^"-2"*")"),
       y = expression("N"["mass_chl"]*" (gN m"^"-2"*")"),
       fill = expression("CO"[2]*" treatment"),
       shape = "Inoculation") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw(base_size = 20)

png("../working_drafts/figs/NxCO2xI_leafN_chl_comps.png",
    width = 10, height = 8, units = 'in', res = 600)
ggarrange(narea.comp, nmass.comp, marea.comp,
          common.legend = TRUE, legend = "right",
          align = "hv", labels = "AUTO")
dev.off()

