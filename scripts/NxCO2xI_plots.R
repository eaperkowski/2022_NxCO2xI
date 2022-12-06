##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(tidyverse)
library(ggpubr)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Load compiled datasheet
df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv",
               na.strings = "NA") %>%
  mutate(n.trt = as.numeric(n.trt)) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nodule.biomass < 0.05)) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE) 

## Add colorblind friendly palette
cbbPalette3 <- c("#DDAA33", "#004488", "#BB5566", "#555555")

##########################################################################
## Ncost regression line prep
##########################################################################
df$ncost[c(101, 102)] <- NA
ncost <- lmer(log(ncost) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(ncost, ~inoc*co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
ncost.full <- data.frame(emmeans(ncost, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

ncost.regline <- data.frame(emmeans(ncost, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(ncost.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"))

##########################################################################
## Ncost plot
##########################################################################
ncost.plot <- ggplot(data = df, aes(x = n.trt, y = ncost)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(ncost.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(ncost.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0,30,10)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
ncost.plot

##########################################################################
## Belowground C regression line prep
##########################################################################
cbg <- lmer(cbg ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
#cbg <- lmer(log(cbg) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

## Emmean fxns for regression lines + error ribbons
cbg.full <- data.frame(emmeans(cbg, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

cbg.regline <- data.frame(emmeans(cbg, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(cbg.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Belowground C plot
##########################################################################
bgc.plot <- ggplot(data = df, aes(x = n.trt, y = cbg)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(cbg.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean), 
              size = 2, se = FALSE) +
  geom_ribbon(data = subset(cbg.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 2, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4,1)) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Belowground carbon biomass (gC)",
       fill = "Treatment", color = "Treatment") +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
bgc.plot

##########################################################################
## Whole plant N regression line prep
##########################################################################
wpn <- lmer(wpn ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

## Emmean fxns for regression lines + error ribbons
wpn.full <- data.frame(emmeans(wpn, ~inoc*co2, "n.trt",
                               at = list(n.trt = seq(0, 630, 5)),
                               type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

wpn.regline <- data.frame(emmeans(wpn, ~1, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(wpn.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Whole plant N plot
##########################################################################
wpn.plot <- ggplot(data = df, aes(x = n.trt, y = wpn)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 4, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(wpn.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean), 
              size = 2, se = FALSE) +
  geom_ribbon(data = subset(wpn.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 2, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 0.48), breaks = seq(0, 0.48, 0.12)) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Whole plant nitrogen biomass (gN)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
wpn.plot

##########################################################################
## Total leaf area regression line prep
##########################################################################
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

## Emmean fxns for regression lines + error ribbons
tla.full <- data.frame(emmeans(tla, ~inoc*co2, "n.trt",
                               at = list(n.trt = seq(0, 630, 5)),
                               type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

tla.regline <- data.frame(emmeans(tla, ~1, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(tla.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Total leaf area plot
##########################################################################
tla.plot <- ggplot(data = df, aes(x = n.trt, y = tla)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(tla.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(tla.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
tla.plot

##########################################################################
## Total biomass regression line prep
##########################################################################
tbio <- lmer(total.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

## Emmean fxns for regression lines + error ribbons
tbio.full <- data.frame(emmeans(tbio, ~inoc*co2, "n.trt",
                               at = list(n.trt = seq(0, 630, 5)),
                               type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

tbio.regline <- data.frame(emmeans(tbio, ~1, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(tbio.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Total biomass plot
##########################################################################
tbio.plot <- ggplot(data = df, aes(x = n.trt, y = total.biomass)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(tbio.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(tbio.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Total biomass (g)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
tbio.plot

##########################################################################
## Root nodule:root biomass regression line prep
##########################################################################
df$nodule.biomass[df$nodule.biomass > 0.05 & df$inoc == "no.inoc"] <- NA
df$nod.root.ratio <- df$nodule.biomass / df$root.biomass

nodroot <- lmer(sqrt(nod.root.ratio) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nodroot, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nodroot.full <- data.frame(emmeans(nodroot, ~inoc*co2, "n.trt",
                                at = list(n.trt = seq(0, 630, 5)),
                                type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

nodroot.regline <- data.frame(emmeans(nodroot, ~1, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(nodroot.full) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc" |
                             co2.inoc == "elv_no.inoc", "dashed", "solid")) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Root nodule:root biomass regression line prep
##########################################################################
nodroot.plot <- ggplot(data = df, aes(x = n.trt, y = nod.root.ratio)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(nodroot.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(nodroot.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Root nodule: root biomass",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  guides(linetype = "none") +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
nodroot.plot

##########################################################################
## Root nodule biomass regression line prep
##########################################################################
df$nodule.biomass[81] <- NA
nod <- lmer(sqrt(nodule.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nod, ~inoc*co2, "n.trt"))


## Emmean fxns for regression lines + error ribbons
nod.full <- data.frame(emmeans(nod, ~inoc*co2, "n.trt",
                               at = list(n.trt = seq(0, 630, 5)),
                               type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

nod.regline <- data.frame(emmeans(nod, ~1, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(nod.full) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc",
                           "dashed", "solid")) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Root nodule biomass regression line prep
##########################################################################
nod.plot <- ggplot(data = df, aes(x = n.trt, y = nodule.biomass)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(nod.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(nod.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.15)) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Root nodule biomass (g)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  guides(linetype = "none") +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
nod.plot

##########################################################################
## Narea regression line prep
##########################################################################
narea <- lmer(narea ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(narea, ~inoc*co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
narea.full <- data.frame(emmeans(narea, ~inoc*co2, "n.trt",
                                at = list(n.trt = seq(0, 630, 5)),
                                type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.regline <- data.frame(emmeans(narea, ~1, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Narea plot
##########################################################################
narea.plot <- ggplot(data = df, aes(x = n.trt, y = narea)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(narea.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(narea.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
narea.plot

##########################################################################
## Nmass regression line prep
##########################################################################
nmass <- lmer(nmass.focal ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nmass, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nmass.full <- data.frame(emmeans(nmass, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

nmass.regline <- data.frame(emmeans(nmass, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(nmass.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"))

##########################################################################
## Nmass plot
##########################################################################
nmass.plot <- ggplot(data = df, aes(x = n.trt, y = nmass.focal)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(nmass.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(nmass.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["mass"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  guides(linetype = "none") +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
nmass.plot

##########################################################################
## Marea regression line prep
##########################################################################
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(marea, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
marea.full <- data.frame(emmeans(marea, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

marea.regline <- data.frame(emmeans(marea, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(marea.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"))

##########################################################################
## Marea plot
##########################################################################
marea.plot <- ggplot(data = df, aes(x = n.trt, y = marea)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(marea.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(marea.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 30)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
marea.plot

##########################################################################
## Chlarea regression line prep
##########################################################################
df$chl.mmolm2[25] <- NA
chlarea <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(chlarea))
car::outlierTest(chlarea)

test(emtrends(chlarea, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
chlarea.full <- data.frame(emmeans(chlarea, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

chlarea.regline <- data.frame(emmeans(chlarea, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(chlarea.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Chl area plot
##########################################################################
chl.plot <- ggplot(data = df, aes(x = n.trt, y = chl.mmolm2)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(chlarea.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(chlarea.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("Chl")["area"]*" (mmol m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
chl.plot

##########################################################################
## Vcmax regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
vcmax25 <- lmer(vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax.full <- data.frame(emmeans(vcmax25, ~inoc*co2, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.regline <- data.frame(emmeans(vcmax25, ~1, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Vcmax plot
##########################################################################
vcmax.plot <- ggplot(data = df, aes(x = n.trt, y = vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(vcmax.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(vcmax.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" (μmol m"^"-2"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
vcmax.plot

##########################################################################
## Jmax regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
jmax25 <- lmer(jmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(jmax25))
test(emtrends(jmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jmax.full <- data.frame(emmeans(jmax25, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

jmax.regline <- data.frame(emmeans(jmax25, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(jmax.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Jmax plot
##########################################################################
jmax.plot <- ggplot(data = df, aes(x = n.trt, y = jmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(jmax.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(jmax.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" (μmol m"^"-2"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
jmax.plot

##########################################################################
## Jmax25:Vcmax25 regression line prep
##########################################################################
df$jmax25.vcmax25[101] <- NA
jvmax <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(jmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jvmax.full <- data.frame(emmeans(jvmax, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

jvmax.regline <- data.frame(emmeans(jvmax, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(jvmax.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Jmax:Vcmax25 plot
##########################################################################
jvmax.plot <- ggplot(data = df, aes(x = n.trt, y = jmax25.vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(jvmax.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(jvmax.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(1.4, 2.2), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*":"*italic(V)["cmax25"])),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
jvmax.plot

##########################################################################
## Rd25 regression line prep
##########################################################################
df$rd25[df$rd25 < 0] <- NA
df$rd25[c(19, 34, 57)] <- NA
rd25 <- lmer(rd25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(rd25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
rd25.full <- data.frame(emmeans(rd25, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

rd25.regline <- data.frame(emmeans(rd25, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(rd25.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Rd25 plot
##########################################################################
rd25.plot <- ggplot(data = df, aes(x = n.trt, y = rd25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(rd25.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(rd25.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("R")["d25"]*" (μmol m"^"-2"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
rd25.plot

##########################################################################
## d25.vcmax25 regression line prep
##########################################################################
df$rd25.vcmax25[df$rd25.vcmax25 < 0] <- NA
df$rd25.vcmax25[c(39, 40, 42)] <- NA
rd25.vcmax25 <- lmer(rd25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

test(emtrends(rd25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
rd25.vcmax25.full <- data.frame(emmeans(rd25.vcmax25, ~inoc*co2, "n.trt",
                                at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

rd25.vcmax25.regline <- data.frame(emmeans(rd25.vcmax25, ~1, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(rd25.vcmax25.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Rd25:Vcmax25 plot
##########################################################################
rd25.vcmax25.plot <- ggplot(data = df, aes(x = n.trt, y = rd25.vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(rd25.vcmax25.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(rd25.vcmax25.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("R")["d25"]*" (μmol m"^"-2"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
rd25.vcmax25.plot

##########################################################################
## Prop leaf N in photosynthesis  regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$p.photo[df$p.photo > 1] <- NA
p.photo <- lmer(p.photo ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(p.photo))
test(emtrends(p.photo, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.photo.full <- data.frame(emmeans(p.photo, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.photo.regline <- data.frame(emmeans(p.photo, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.photo.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc", 
                           "dashed", "solid"))

##########################################################################
## Prop leaf N in photosynthesis plot
##########################################################################
p.photo.plot <- ggplot(data = df, aes(x = n.trt, y = p.photo)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(p.photo.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(p.photo.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["photo"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
p.photo.plot

##########################################################################
## Prop leaf N in rubisco regression line prep
##########################################################################
df$p.rubisco[c(45)] <- NA
p.rub <- lmer(p.rubisco ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(p.rub))
test(emtrends(p.rub, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.rub.full <- data.frame(emmeans(p.rub, ~inoc*co2, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.rub.regline <- data.frame(emmeans(p.rub, ~1, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.rub.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_no.inoc", 
                           "dashed", "solid"))

##########################################################################
## Prop leaf N in rubisco plot
##########################################################################
p.rub.plot <- ggplot(data = df, aes(x = n.trt, y = p.rubisco)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(p.rub.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(p.rub.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["rubisco"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
p.rub.plot

##########################################################################
## Prop leaf N in bioenergetics regression line prep
##########################################################################
df$p.bioe[c(45)] <- NA
p.bioe <- lmer(p.bioe ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(p.bioe))
test(emtrends(p.bioe, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.bioe.full <- data.frame(emmeans(p.bioe, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)))) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.bioe.regline <- data.frame(emmeans(p.bioe, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)))) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.bioe.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_no.inoc", 
                           "dashed", "solid"))

##########################################################################
## Prop leaf N in bioenergetics plot
##########################################################################
p.bioe.plot <- ggplot(data = df, aes(x = n.trt, y = p.bioe)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(p.bioe.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(p.bioe.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.03)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["bioenergetics"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
p.bioe.plot

##########################################################################
## Prop leaf N in light harvesting proteins regression line prep
##########################################################################
df$p.lightharv[c(25, 39, 45)] <- NA
p.light <- lmer(p.lightharv ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(p.light))
test(emtrends(p.light, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.light.full <- data.frame(emmeans(p.light, ~inoc*co2, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.light.regline <- data.frame(emmeans(p.light, ~1, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.light.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc", 
                           "dashed", "solid"))

##########################################################################
## Prop leaf N in light harvesting proteins plot
##########################################################################
p.light.plot <- ggplot(data = df, aes(x = n.trt, y = p.lightharv)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(p.light.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(p.light.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["light"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
p.light.plot

##########################################################################
## Prop leaf N in structure regression line prep
##########################################################################
df$p.structure[c(39, 45, 101, 102, 104)] <- NA
p.str <- lmer(log(p.structure) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(p.str))
test(emtrends(p.str, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.str.full <- data.frame(emmeans(p.str, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.str.regline <- data.frame(emmeans(p.str, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.str.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_no.inoc", 
                           "dashed", "solid")) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc", 
                           "dashed", "solid"))

##########################################################################
## Prop leaf N in structure plot
##########################################################################
p.str.plot <- ggplot(data = df, aes(x = n.trt, y = p.structure)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(p.str.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(p.str.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["structure"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
p.str.plot

##########################################################################
## PNUE regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$pnue[45] <- NA
pnue <- lmer(pnue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(pnue))
test(emtrends(pnue, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
pnue.full <- data.frame(emmeans(pnue, ~inoc*co2, "n.trt",
                                 at = list(n.trt = seq(0, 630, 5)),
                                 type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

pnue.regline <- data.frame(emmeans(pnue, ~1, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(pnue.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_no.inoc", 
                           "dashed", "solid"))

##########################################################################
## PNUE plot
##########################################################################
pnue.plot <- ggplot(data = df, aes(x = n.trt, y = pnue)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(pnue.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(pnue.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("PNUE (μmol CO"["2"]*" gN"^"-1"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
pnue.plot

##########################################################################
## iWUE regression line prep - to be replaced with chi?
##########################################################################
iwue <- lmer(iwue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(iwue))
test(emtrends(iwue, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
iwue.full <- data.frame(emmeans(iwue, ~inoc*co2, "n.trt",
                                at = list(n.trt = seq(0, 630, 5)),
                                type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

iwue.regline <- data.frame(emmeans(iwue, ~1, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(iwue.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## iWUE plot
##########################################################################
iwue.plot <- ggplot(data = df, aes(x = n.trt, y = iwue)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(iwue.regline, co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(iwue.regline, co2.inoc == "overall"),
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(25, 125), breaks = seq(25, 125, 25)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("iWUE (μmol CO"["2"]*" mol"^"-1"*"H"["2"]*"O)")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14))
iwue.plot

##########################################################################
## Narea:gs regression line prep
##########################################################################
narea.gs <- lmer(log(narea.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(narea.gs))
test(emtrends(narea.gs, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
narea.gs.full <- data.frame(emmeans(narea.gs, ~inoc*co2, "n.trt",
                                at = list(n.trt = seq(0, 630, 5)),
                                type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.gs.regline <- data.frame(emmeans(narea.gs, ~1, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.gs.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc" | co2.inoc == "elv_no.inoc", 
                           "dashed", "solid"))

##########################################################################
## Narea:gs plot
##########################################################################
narea.gs.plot <- ggplot(data = df, aes(x = n.trt, y = narea.gs)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(narea.gs.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(narea.gs.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 32, 8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*": "*italic("g")["s"]*" (gN s mol"^"-1"*" H"["2"]*"O)")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
narea.gs.plot

##########################################################################
## Vcmax:gs regression line prep
##########################################################################
vcmax.gs <- lmer(log(vcmax.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(vcmax.gs))
test(emtrends(vcmax.gs, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax.gs.full <- data.frame(emmeans(vcmax.gs, ~inoc*co2, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.gs.regline <- data.frame(emmeans(vcmax.gs, ~1, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.gs.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Vcmax:gs plot
##########################################################################
vcmax.gs.plot <- ggplot(data = df, aes(x = n.trt, y = vcmax.gs)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(vcmax.gs.regline, co2.inoc == "overall"),
              aes(y = response), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(vcmax.gs.regline, co2.inoc == "overall"),
              aes(y = response, ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 200)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*": "*italic("g")["s"]*" (μmol CO"["2"]*" mol"^"-1"*" H"["2"]*"O)")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14))
vcmax.gs.plot

##########################################################################
## Figure 1: whole plant plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig1_wholePlant.png",
    height = 12, width = 7.5, units = "in", res = 600)
ggarrange(ncost.plot, tla.plot, tbio.plot, ncol = 1, nrow = 3,
          common.legend = TRUE, align = "hv",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 2: nitrogen fixation
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig2_nfix.png",
    height = 4.5, width = 12, units = "in", res = 600)
ggarrange(nod.plot, nodroot.plot, ncol = 2, nrow = 1,
          common.legend = TRUE, align = "hv",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 3: leaf N plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig3_leafN.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(narea.plot, nmass.plot, marea.plot, chl.plot,
          ncol = 2, nrow = 2,
          common.legend = TRUE, align = "hv",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 4: leaf physiology plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig4_photo.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(vcmax.plot, jmax.plot, jvmax.plot, rd25.plot, ncol = 2, nrow = 2,
          common.legend = TRUE, align = "hv",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 5: propN photosynthesis/structure
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig5_propN.png",
    height = 8.5, width = 17, units = "in", res = 600)
ggarrange(p.rub.plot, p.bioe.plot, p.light.plot,
          p.photo.plot,
          p.str.plot, ncol = 3, nrow = 2,
          common.legend = TRUE, align = "hv",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 5: PNUE/iWUE tradeoffs
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig6_PNUE_iWUE.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(pnue.plot, iwue.plot, narea.gs.plot, vcmax.gs.plot,
          ncol = 2, nrow = 2, common.legend = TRUE, align = "hv",
          legend = "right", labels = "AUTO",
          font.label = list(size = 18))
dev.off()







