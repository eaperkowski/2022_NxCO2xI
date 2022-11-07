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
  

## Do I need to roll all chlorophyll leaves?
narea.comp <- ggplot(data = df, aes(x = narea, y = narea.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc)) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_cor(label.y = 2.9) +
  stat_regline_equation(label.y = 2.8) +
  scale_shape_manual(values = c(21, 22)) +
  scale_x_continuous(limits = c(0,3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(0,3), breaks = seq(0, 3, 1)) +
  labs(x = expression("N"["area"]*" (gN m"^"-2"*")"),
       y = expression("N"["area_chl"]*" (gN m"^"-2"*")")) +
  theme_bw(base_size = 20)

png("../working_drafts/figs/NxCO2xI_narea_chl_comp.png",
    width = 8, height = 8, units = 'in', res = 600)
narea.comp
dev.off()

marea.comp <- ggplot(data = subset(df, marea.chl < 100), 
       aes(x = marea, y = marea.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc)) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_cor(label.y = 100) +
  stat_regline_equation(label.y = 97) +
  scale_shape_manual(values = c(21, 22)) +
  scale_x_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  scale_y_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  labs(x = expression("M"["area"]*" (g m"^"-2"*")"),
       y = expression("M"["area_chl"]*" (g m"^"-2"*")")) +
  theme_bw(base_size = 18)

png("../working_drafts/figs/NxCO2xI_marea_chl_comp.png",
    width = 8, height = 8, units = 'in', res = 600)
marea.comp
dev.off()


nmass.comp <- ggplot(data = df, aes(x = nmass.focal, y = nmass.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc)) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_cor(label.y = 0.075) +
  stat_regline_equation(label.y = 0.069) +
  scale_shape_manual(values = c(21, 22)) +
  scale_x_continuous(limits = c(0, 0.075), breaks = seq(0, 0.075, 0.025)) +
  scale_y_continuous(limits = c(0, 0.075), breaks = seq(0, 0.075, 0.025)) +
  labs(x = expression("N"["mass"]*" (gN m"^"-2"*")"),
       y = expression("N"["mass_chl"]*" (gN m"^"-2"*")")) +
  theme_bw(base_size = 20)

png("../working_drafts/NxCO2xI_nmass_chl_comp.png",
    width = 8, height = 8, units = 'in', res = 600)
nmass.comp
dev.off()

## Preliminary models for LEMONTREE meeting
vcmax25 <- lm(vcmax25 ~ as.factor(co2) * inoc * as.numeric(n.trt),
              data = df)
summary(vcmax25)
Anova(vcmax25)

shapiro.test(residuals(vcmax25))
outlierTest(vcmax25)

## Preliminary plots for LEMONTREE meeting
inoc.labs <- c("Not inoculated", "Inoculated")
names(inoc.labs) <- c("n", "y")


vcmax25 <- ggplot(data = subset(df, !is.na(inoc)), 
       aes(x = as.numeric(n.trt), y = vcmax25)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 30)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  labs(x = NULL,
       y = expression(italic("V")["cmax25"]*" (μmol m"^"-2"*"s"^"-1"*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 22)

jmax25 <- ggplot(data = subset(df, !is.na(inoc)), 
       aes(x = as.numeric(n.trt), y = jmax25)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("J")["max25"]*" (μmol m"^"-2"*"s"^"-1"*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  theme_bw(base_size = 22)


jv.rat <- ggplot(data = subset(df, !is.na(inoc)), 
       aes(x = as.numeric(n.trt), y = jmax25.vcmax25)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  scale_y_continuous(limits = c(1.4,2), breaks = seq(1.4,2,0.2)) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("J")["max25"]*":"*italic("V")["cmax25"]),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 22)

tla <- ggplot(data = subset(df, !is.na(inoc)), 
              aes(x = as.numeric(n.trt), y = tla)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 300)) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression("Total leaf area (cm"^2*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 22)

png("../working_drafts/NxCO2_LEMONTREE_fig_vcmax_jvrat.png",
    width = 12, height = 12, units = 'in', res = 600)
ggarrange(vcmax25, jv.rat, nrow = 2, align = "hv", common.legend = TRUE,
          legend = "right")
dev.off()

png("../working_drafts/NxCO2_LEMONTREE_fig_tla.png",
    width = 12, height = 6, units = 'in', res = 600)
tla
dev.off()

