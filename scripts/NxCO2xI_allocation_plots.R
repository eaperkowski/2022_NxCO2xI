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

## Calculate RMF, SMF, LMF
df$rmf <- df$root.biomass / df$total.biomass
df$smf <- df$stem.biomass / df$total.biomass
df$lmf <- df$leaf.biomass / df$total.biomass

##########################################################################
## Root mass fraction
##########################################################################
df$rmf[102] <- NA

## model
rmf <- lmer(log(rmf) ~ n.trt * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(rmf))

## trends
test(emtrends(rmf, ~co2*inoc, "n.trt"))
## no n.trt effect in inoculated pots, rmf reduction in uninoculated pots

## Emmean fxns for regression lines + error ribbons
rmf.full <- data.frame(emmeans(rmf, ~inoc*co2, "n.trt",
                                       at = list(n.trt = seq(0, 630, 1)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

rmf.regline <- data.frame(emmeans(rmf, ~1, "n.trt",
                                          at = list(n.trt = seq(0, 630, 1)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(rmf.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

## Plot
rmf.plot <- ggplot(data = df, aes(x = n.trt, y = rmf)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(rmf.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(rmf.regline, co2.inoc != "overall"),
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
  scale_y_continuous(limits = c(0.2, 0.62), breaks = seq(0.2, 0.6, 0.1)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Root mass fraction (g"["root"]*" g"["total"]*""^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
rmf.plot

##########################################################################
## Stem mass fraction
##########################################################################
df$smf[114] <- NA

## model
smf <- lmer(smf ~ n.trt * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(smf))

## trends
test(emtrends(smf, ~co2*inoc, "n.trt"))
## no n.trt effect in inoculated pots, rmf reduction in uninoculated pots

## Emmean fxns for regression lines + error ribbons
smf.full <- data.frame(emmeans(smf, ~inoc*co2, "n.trt",
                               at = list(n.trt = seq(0, 630, 1)),
                               type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

smf.regline <- data.frame(emmeans(smf, ~1, "n.trt",
                                  at = list(n.trt = seq(0, 630, 1)),
                                  type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(smf.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "solid", "dashed"))

## Plot
smf.plot <- ggplot(data = df, aes(x = n.trt, y = smf)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(smf.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(smf.regline, co2.inoc != "overall"),
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
  scale_y_continuous(limits = c(0.1, 0.3), breaks = seq(0.1, 0.3, 0.05)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Stem mass fraction (g"["stem"]*" g"["total"]*""^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
smf.plot

##########################################################################
## Leaf mass fraction
##########################################################################
## model
lmf <- lmer(lmf ~ n.trt * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(lmf))
car::outlierTest(lmf)

## trends
test(emtrends(lmf, ~co2*inoc, "n.trt"))
## no n.trt effect in inoculated pots, rmf reduction in uninoculated pots

## Emmean fxns for regression lines + error ribbons
lmf.full <- data.frame(emmeans(lmf, ~inoc*co2, "n.trt",
                               at = list(n.trt = seq(0, 630, 1)),
                               type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

lmf.regline <- data.frame(emmeans(lmf, ~1, "n.trt",
                                  at = list(n.trt = seq(0, 630, 1)),
                                  type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(lmf.full) %>%
  dplyr::select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "solid", "dashed"))

## Plot
lmf.plot <- ggplot(data = df, aes(x = n.trt, y = lmf)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(lmf.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(lmf.regline, co2.inoc != "overall"),
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
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_y_continuous(limits = c(0.1, 0.51), breaks = seq(0.1, 0.5, 0.1)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Leaf mass fraction (g"["leaf"]*" g"["total"]*""^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
lmf.plot

