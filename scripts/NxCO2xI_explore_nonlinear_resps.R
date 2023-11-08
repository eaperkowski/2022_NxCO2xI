##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)
library(nlme)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Read in compiled data file
df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  mutate(n.trt = as.numeric(n.trt),
         inoc = factor(inoc, levels = c("no.inoc", "inoc")),
         co2 = factor(co2, levels = c("amb", "elv")),
         co2.inoc = str_c(co2, "_", inoc),
         nod.root.ratio = nodule.biomass / root.biomass,
         pnue.growth = anet.growth / (narea / 14)) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05))
## filter all uninoculated pots that have nod biomass > 0.05 g;
## hard code inoc/co2 to make coefficients easier to understand

df.removed <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  mutate(nod.root.ratio = nodule.biomass / root.biomass) %>%
  dplyr::filter(inoc == "no.inoc" & nod.root.ratio > 0.05)

##########################################################################
## Vcmax25 nonlinear saturating regression
##########################################################################

## Create nonlinear saturating regressions
vcmax.nls.elv <- nls(formula = vcmax25 ~ a + ((b*n.trt)/(c+n.trt)),
                 data = subset(df, inoc == "no.inoc" & co2 == "elv"),
                 start = list(a = 18, b = 0.3, c = -0.001))

vcmax.nls.amb <- nls(formula = vcmax25 ~ a + ((b*n.trt)/(c+n.trt)),
                     data = subset(df, inoc == "no.inoc" & co2 == "amb"),
                     start = list(a = 18, b = 0.3, c = -0.001))

## Create predicted trendlines for both CO2 treatments
vcmax.nls.elv.pred <- data.frame(emmeans(vcmax.nls.elv, ~1, "n.trt", 
                                     at = list(n.trt = seq(0, 630, 1)),
                                     data = subset(df, inoc == "no.inoc")))
vcmax.nls.amb.pred <- data.frame(emmeans(vcmax.nls.amb, ~1, "n.trt", 
                                         at = list(n.trt = seq(0, 630, 1)),
                                         data = subset(df, inoc == "no.inoc")))

## Create Vcmax plot
vcmax25.plot <- ggplot(data = subset(df, inoc == "no.inoc"), 
                       aes(x = n.trt, y = vcmax25)) +
  geom_point(aes(fill = co2), alpha = 0.75, size = 4, shape = 21) +
  geom_smooth(data = vcmax.nls.elv.pred, aes(y = emmean), color = "#b2182b",
              linewidth = 2, se = FALSE) +
  geom_smooth(data = vcmax.nls.amb.pred, aes(y = emmean), color = "#2166ac",
              linewidth = 2, se = FALSE) +
  scale_fill_manual(values = c("#2166ac", "#b2182b"),
                     labels = c(expression("Ambient CO"["2"]),
                                expression("Elevated CO"["2"]))) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 30)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bold("CO"["2"]*" treatment"))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
vcmax25.plot

##########################################################################
## Jmax25 nonlinear saturating regression
##########################################################################

## Create nonlinear saturating regressions
jmax.nls.elv <- nls(formula = jmax25 ~ a + ((b*n.trt)/(c + n.trt)),
                     data = subset(df, inoc == "no.inoc" & co2 == "elv"),
                     start = list(a = 46, b = 185, c = 347))

jmax.nls.amb <- nls(formula = jmax25 ~ a + ((b*n.trt)/(c+n.trt)),
                     data = subset(df, inoc == "no.inoc" & co2 == "amb"),
                     start = list(a = 46, b = 185, c = 347))

## Create predicted trendlines for both CO2 treatments
jmax.nls.elv.pred <- data.frame(emmeans(jmax.nls.elv, ~1, "n.trt", 
                                         at = list(n.trt = seq(0, 630, 1)),
                                         data = subset(df, inoc == "no.inoc")))
jmax.nls.amb.pred <- data.frame(emmeans(jmax.nls.amb, ~1, "n.trt", 
                                         at = list(n.trt = seq(0, 630, 1)),
                                         data = subset(df, inoc == "no.inoc")))

## Create Jmax plot
jmax25.plot <- ggplot(data = subset(df, inoc == "no.inoc"), 
                       aes(x = n.trt, y = jmax25)) +
  geom_point(aes(fill = co2), alpha = 0.75, size = 4, shape = 21) +
  geom_smooth(data = jmax.nls.elv.pred, aes(y = emmean), color = "#b2182b",
              linewidth = 2, se = FALSE) +
  geom_smooth(data = jmax.nls.amb.pred, aes(y = emmean), color = "#2166ac",
              linewidth = 2, se = FALSE) +
  scale_fill_manual(values = c("#2166ac", "#b2182b"),
                    labels = c(expression("Ambient CO"["2"]),
                               expression("Elevated CO"["2"]))) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bold("CO"["2"]*" treatment"))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
jmax25.plot

##########################################################################
## Write merged Vcmax25 and Jmax25 plot
##########################################################################
png("../working_drafts/figs/NxCO2xI_photo_nonlinear.png", width = 12, 
    height = 4.5, units = "in", res = 600)
ggarrange(vcmax25.plot, jmax25.plot, common.legend = TRUE, legend = "right",
          labels = c("a", "b"))
dev.off()


