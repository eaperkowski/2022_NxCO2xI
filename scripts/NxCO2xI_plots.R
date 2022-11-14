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
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

## Add colorblind friendly palette
cbbPalette3 <- c("#DDAA33", "#004488", "#BB5566", "#555555")

##########################################################################
## Ncost regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$ncost[c(110,111)] <- NA
ncost <- lmer(log(ncost) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(ncost))
Anova(ncost)

## Emmean fxns for regression lines + error ribbons
ncost.full.0 <- data.frame(n.trt = 0, emmeans(ncost, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
ncost.full.70 <- data.frame(n.trt = 70, emmeans(ncost, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
ncost.full.140 <- data.frame(n.trt = 140, emmeans(ncost, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
ncost.full.210 <- data.frame(n.trt = 210, emmeans(ncost, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
ncost.full.280 <- data.frame(n.trt = 280, emmeans(ncost, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
ncost.full.350 <- data.frame(n.trt = 350, emmeans(ncost, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
ncost.full.630 <- data.frame(n.trt = 630, emmeans(ncost, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
ncost.full <- ncost.full.0 %>% full_join(ncost.full.70) %>% 
  full_join(ncost.full.140) %>% full_join(ncost.full.210) %>% 
  full_join(ncost.full.280) %>% full_join(ncost.full.350) %>% 
  full_join(ncost.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

ncost.ntrt.0 <- data.frame(n.trt = 0, emmeans(ncost, ~1, 
                                              at = list(n.trt = 0), type = "response"))
ncost.ntrt.70 <- data.frame(n.trt = 70, emmeans(ncost, ~1, 
                                               at = list(n.trt = 70), type = "response"))
ncost.ntrt.140 <- data.frame(n.trt = 140, emmeans(ncost, ~1, 
                                                at = list(n.trt = 140), type = "response"))
ncost.ntrt.210 <- data.frame(n.trt = 210, emmeans(ncost, ~1, 
                                                at = list(n.trt = 210), type = "response"))
ncost.ntrt.280 <- data.frame(n.trt = 280, emmeans(ncost, ~1, 
                                                at = list(n.trt = 280), type = "response"))
ncost.ntrt.350 <- data.frame(n.trt = 350, emmeans(ncost, ~1, 
                                                at = list(n.trt = 350), type = "response"))
ncost.ntrt.630 <- data.frame(n.trt = 630, emmeans(ncost, ~1, 
                                                at = list(n.trt = 630), type = "response"))
ncost.regline <- ncost.ntrt.0 %>% full_join(ncost.ntrt.70) %>%
  full_join(ncost.ntrt.140) %>% full_join(ncost.ntrt.210) %>%
  full_join(ncost.ntrt.280) %>% full_join(ncost.ntrt.350) %>%
  full_join(ncost.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(ncost.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Copy removed outliers and lmer fxn
df$cbg[c(68, 70, 61)] <- NA

cbg <- lmer(cbg ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(cbg))
outlierTest(cbg)

## Emmean fxns for regression lines + error ribbons
cbg.full.0 <- data.frame(n.trt = 0, emmeans(cbg, ~inoc*co2, 
                                            at = list(n.trt = 0), 
                                            type = "response"))
cbg.full.70 <- data.frame(n.trt = 70, emmeans(cbg, ~inoc*co2, 
                                              at = list(n.trt = 70), 
                                              type = "response"))
cbg.full.140 <- data.frame(n.trt = 140, emmeans(cbg, ~inoc*co2, 
                                                at = list(n.trt = 140), 
                                                type = "response")) 
cbg.full.210 <- data.frame(n.trt = 210, emmeans(cbg, ~inoc*co2, 
                                                at = list(n.trt = 210), 
                                                type = "response"))
cbg.full.280 <- data.frame(n.trt = 280, emmeans(cbg, ~inoc*co2, 
                                                at = list(n.trt = 280), 
                                                type = "response")) 
cbg.full.350 <- data.frame(n.trt = 350, emmeans(cbg, ~inoc*co2, 
                                                at = list(n.trt = 350), 
                                                type = "response")) 
cbg.full.630 <- data.frame(n.trt = 630, emmeans(cbg, ~inoc*co2, 
                                                at = list(n.trt = 630),
                                                type = "response")) 
cbg.full <- cbg.full.0 %>% full_join(cbg.full.70) %>% 
  full_join(cbg.full.140) %>% full_join(cbg.full.210) %>% 
  full_join(cbg.full.280) %>% full_join(cbg.full.350) %>% 
  full_join(cbg.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

cbg.ntrt.0 <- data.frame(n.trt = 0, emmeans(cbg, ~1, 
                                            at = list(n.trt = 0), type = "response"))
cbg.ntrt.70 <- data.frame(n.trt = 70, emmeans(cbg, ~1, 
                                              at = list(n.trt = 70), type = "response"))
cbg.ntrt.140 <- data.frame(n.trt = 140, emmeans(cbg, ~1, 
                                                at = list(n.trt = 140), type = "response"))
cbg.ntrt.210 <- data.frame(n.trt = 210, emmeans(cbg, ~1, 
                                                at = list(n.trt = 210), type = "response"))
cbg.ntrt.280 <- data.frame(n.trt = 280, emmeans(cbg, ~1, 
                                                at = list(n.trt = 280), type = "response"))
cbg.ntrt.350 <- data.frame(n.trt = 350, emmeans(cbg, ~1, 
                                                at = list(n.trt = 350), type = "response"))
cbg.ntrt.630 <- data.frame(n.trt = 630, emmeans(cbg, ~1, 
                                                at = list(n.trt = 630), type = "response"))
cbg.regline <- cbg.ntrt.0 %>% full_join(cbg.ntrt.70) %>%
  full_join(cbg.ntrt.140) %>% full_join(cbg.ntrt.210) %>%
  full_join(cbg.ntrt.280) %>% full_join(cbg.ntrt.350) %>%
  full_join(cbg.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(cbg.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4,1)) +
  labs(x = "Soil N fertilization (ppm twice per week)",
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
## Copy removed outliers and lmer fxn
wpn <- lmer(sqrt(wpn) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(wpn)
test(emtrends(wpn, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
wpn.full.0 <- data.frame(n.trt = 0, emmeans(wpn, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
wpn.full.70 <- data.frame(n.trt = 70, emmeans(wpn, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
wpn.full.140 <- data.frame(n.trt = 140, emmeans(wpn, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
wpn.full.210 <- data.frame(n.trt = 210, emmeans(wpn, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
wpn.full.280 <- data.frame(n.trt = 280, emmeans(wpn, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
wpn.full.350 <- data.frame(n.trt = 350, emmeans(wpn, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
wpn.full.630 <- data.frame(n.trt = 630, emmeans(wpn, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
wpn.full <- wpn.full.0 %>% full_join(wpn.full.70) %>% 
  full_join(wpn.full.140) %>% full_join(wpn.full.210) %>% 
  full_join(wpn.full.280) %>% full_join(wpn.full.350) %>% 
  full_join(wpn.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

wpn.ntrt.0 <- data.frame(n.trt = 0, emmeans(wpn, ~1, 
                                              at = list(n.trt = 0), type = "response"))
wpn.ntrt.70 <- data.frame(n.trt = 70, emmeans(wpn, ~1, 
                                                at = list(n.trt = 70), type = "response"))
wpn.ntrt.140 <- data.frame(n.trt = 140, emmeans(wpn, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
wpn.ntrt.210 <- data.frame(n.trt = 210, emmeans(wpn, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
wpn.ntrt.280 <- data.frame(n.trt = 280, emmeans(wpn, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
wpn.ntrt.350 <- data.frame(n.trt = 350, emmeans(wpn, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
wpn.ntrt.630 <- data.frame(n.trt = 630, emmeans(wpn, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
wpn.regline <- wpn.ntrt.0 %>% full_join(wpn.ntrt.70) %>%
  full_join(wpn.ntrt.140) %>% full_join(wpn.ntrt.210) %>%
  full_join(wpn.ntrt.280) %>% full_join(wpn.ntrt.350) %>%
  full_join(wpn.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(wpn.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Whole plant N plot
##########################################################################
wpn.plot <- ggplot(data = df, aes(x = n.trt, y = wpn)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 4, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(wpn.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response), 
              size = 2, se = FALSE) +
  geom_ribbon(data = subset(wpn.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 2, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 0.48), breaks = seq(0, 0.48, 0.12)) +
  labs(x = "Soil N fertilization (ppm twice per week)",
       y = expression(bold("Carbon cost to acquire nitrogen (gC gN"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
wpn.plot

##########################################################################
## Total leaf area regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(tla)
shapiro.test(residuals(tla))
test(emtrends(tla, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
tla.full.0 <- data.frame(n.trt = 0, emmeans(tla, ~inoc*co2, 
                                            at = list(n.trt = 0), 
                                            type = "response"))
tla.full.70 <- data.frame(n.trt = 70, emmeans(tla, ~inoc*co2, 
                                              at = list(n.trt = 70), 
                                              type = "response"))
tla.full.140 <- data.frame(n.trt = 140, emmeans(tla, ~inoc*co2, 
                                                at = list(n.trt = 140), 
                                                type = "response")) 
tla.full.210 <- data.frame(n.trt = 210, emmeans(tla, ~inoc*co2, 
                                                at = list(n.trt = 210), 
                                                type = "response"))
tla.full.280 <- data.frame(n.trt = 280, emmeans(tla, ~inoc*co2, 
                                                at = list(n.trt = 280), 
                                                type = "response")) 
tla.full.350 <- data.frame(n.trt = 350, emmeans(tla, ~inoc*co2, 
                                                at = list(n.trt = 350), 
                                                type = "response")) 
tla.full.630 <- data.frame(n.trt = 630, emmeans(tla, ~inoc*co2, 
                                                at = list(n.trt = 630),
                                                type = "response")) 
tla.full <- tla.full.0 %>% full_join(tla.full.70) %>% 
  full_join(tla.full.140) %>% full_join(tla.full.210) %>% 
  full_join(tla.full.280) %>% full_join(tla.full.350) %>% 
  full_join(tla.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

tla.ntrt.0 <- data.frame(n.trt = 0, emmeans(tla, ~1, 
                                            at = list(n.trt = 0), type = "response"))
tla.ntrt.70 <- data.frame(n.trt = 70, emmeans(tla, ~1, 
                                              at = list(n.trt = 70), type = "response"))
tla.ntrt.140 <- data.frame(n.trt = 140, emmeans(tla, ~1, 
                                                at = list(n.trt = 140), type = "response"))
tla.ntrt.210 <- data.frame(n.trt = 210, emmeans(tla, ~1, 
                                                at = list(n.trt = 210), type = "response"))
tla.ntrt.280 <- data.frame(n.trt = 280, emmeans(tla, ~1, 
                                                at = list(n.trt = 280), type = "response"))
tla.ntrt.350 <- data.frame(n.trt = 350, emmeans(tla, ~1, 
                                                at = list(n.trt = 350), type = "response"))
tla.ntrt.630 <- data.frame(n.trt = 630, emmeans(tla, ~1, 
                                                at = list(n.trt = 630), type = "response"))
tla.regline <- tla.ntrt.0 %>% full_join(tla.ntrt.70) %>%
  full_join(tla.ntrt.140) %>% full_join(tla.ntrt.210) %>%
  full_join(tla.ntrt.280) %>% full_join(tla.ntrt.350) %>%
  full_join(tla.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(tla.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))

##########################################################################
## Total biomass regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
tbio <- lmer(total.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(tbio)
shapiro.test(residuals(tbio))
outlierTest(tbio)

test(emtrends(tbio, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
tbio.full.0 <- data.frame(n.trt = 0, emmeans(tbio, ~inoc*co2, 
                                            at = list(n.trt = 0), 
                                            type = "response"))
tbio.full.70 <- data.frame(n.trt = 70, emmeans(tbio, ~inoc*co2, 
                                              at = list(n.trt = 70), 
                                              type = "response"))
tbio.full.140 <- data.frame(n.trt = 140, emmeans(tbio, ~inoc*co2, 
                                                at = list(n.trt = 140), 
                                                type = "response")) 
tbio.full.210 <- data.frame(n.trt = 210, emmeans(tbio, ~inoc*co2, 
                                                at = list(n.trt = 210), 
                                                type = "response"))
tbio.full.280 <- data.frame(n.trt = 280, emmeans(tbio, ~inoc*co2, 
                                                at = list(n.trt = 280), 
                                                type = "response")) 
tbio.full.350 <- data.frame(n.trt = 350, emmeans(tbio, ~inoc*co2, 
                                                at = list(n.trt = 350), 
                                                type = "response")) 
tbio.full.630 <- data.frame(n.trt = 630, emmeans(tbio, ~inoc*co2, 
                                                at = list(n.trt = 630),
                                                type = "response")) 
tbio.full <- tbio.full.0 %>% full_join(tbio.full.70) %>% 
  full_join(tbio.full.140) %>% full_join(tbio.full.210) %>% 
  full_join(tbio.full.280) %>% full_join(tbio.full.350) %>% 
  full_join(tbio.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

tbio.ntrt.0 <- data.frame(n.trt = 0, emmeans(tbio, ~1, 
                                            at = list(n.trt = 0), type = "response"))
tbio.ntrt.70 <- data.frame(n.trt = 70, emmeans(tbio, ~1, 
                                              at = list(n.trt = 70), type = "response"))
tbio.ntrt.140 <- data.frame(n.trt = 140, emmeans(tbio, ~1, 
                                                at = list(n.trt = 140), type = "response"))
tbio.ntrt.210 <- data.frame(n.trt = 210, emmeans(tbio, ~1, 
                                                at = list(n.trt = 210), type = "response"))
tbio.ntrt.280 <- data.frame(n.trt = 280, emmeans(tbio, ~1, 
                                                at = list(n.trt = 280), type = "response"))
tbio.ntrt.350 <- data.frame(n.trt = 350, emmeans(tbio, ~1, 
                                                at = list(n.trt = 350), type = "response"))
tbio.ntrt.630 <- data.frame(n.trt = 630, emmeans(tbio, ~1, 
                                                at = list(n.trt = 630), type = "response"))
tbio.regline <- tbio.ntrt.0 %>% full_join(tbio.ntrt.70) %>%
  full_join(tbio.ntrt.140) %>% full_join(tbio.ntrt.210) %>%
  full_join(tbio.ntrt.280) %>% full_join(tbio.ntrt.350) %>%
  full_join(tbio.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(tbio.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Narea regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
narea <- lmer(narea ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(narea)
shapiro.test(residuals(narea))
outlierTest(narea)

test(emtrends(narea, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
narea.full.0 <- data.frame(n.trt = 0, emmeans(narea, ~inoc*co2, 
                                             at = list(n.trt = 0), 
                                             type = "response"))
narea.full.70 <- data.frame(n.trt = 70, emmeans(narea, ~inoc*co2, 
                                               at = list(n.trt = 70), 
                                               type = "response"))
narea.full.140 <- data.frame(n.trt = 140, emmeans(narea, ~inoc*co2, 
                                                 at = list(n.trt = 140), 
                                                 type = "response")) 
narea.full.210 <- data.frame(n.trt = 210, emmeans(narea, ~inoc*co2, 
                                                 at = list(n.trt = 210), 
                                                 type = "response"))
narea.full.280 <- data.frame(n.trt = 280, emmeans(narea, ~inoc*co2, 
                                                 at = list(n.trt = 280), 
                                                 type = "response")) 
narea.full.350 <- data.frame(n.trt = 350, emmeans(narea, ~inoc*co2, 
                                                 at = list(n.trt = 350), 
                                                 type = "response")) 
narea.full.630 <- data.frame(n.trt = 630, emmeans(narea, ~inoc*co2, 
                                                 at = list(n.trt = 630),
                                                 type = "response")) 
narea.full <- narea.full.0 %>% full_join(narea.full.70) %>% 
  full_join(narea.full.140) %>% full_join(narea.full.210) %>% 
  full_join(narea.full.280) %>% full_join(narea.full.350) %>% 
  full_join(narea.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.ntrt.0 <- data.frame(n.trt = 0, emmeans(narea, ~1, 
                                             at = list(n.trt = 0), type = "response"))
narea.ntrt.70 <- data.frame(n.trt = 70, emmeans(narea, ~1, 
                                               at = list(n.trt = 70), type = "response"))
narea.ntrt.140 <- data.frame(n.trt = 140, emmeans(narea, ~1, 
                                                 at = list(n.trt = 140), type = "response"))
narea.ntrt.210 <- data.frame(n.trt = 210, emmeans(narea, ~1, 
                                                 at = list(n.trt = 210), type = "response"))
narea.ntrt.280 <- data.frame(n.trt = 280, emmeans(narea, ~1, 
                                                 at = list(n.trt = 280), type = "response"))
narea.ntrt.350 <- data.frame(n.trt = 350, emmeans(narea, ~1, 
                                                 at = list(n.trt = 350), type = "response"))
narea.ntrt.630 <- data.frame(n.trt = 630, emmeans(narea, ~1, 
                                                 at = list(n.trt = 630), type = "response"))
narea.regline <- narea.ntrt.0 %>% full_join(narea.ntrt.70) %>%
  full_join(narea.ntrt.140) %>% full_join(narea.ntrt.210) %>%
  full_join(narea.ntrt.280) %>% full_join(narea.ntrt.350) %>%
  full_join(narea.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Copy removed outliers and lmer fxn
df$nmass.focal[c(39, 50, 110, 111, 114)] <- NA
nmass <- lmer(log(nmass.focal) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(nmass)
shapiro.test(residuals(nmass))
outlierTest(nmass)

test(emtrends(nmass, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
nmass.full.0 <- data.frame(n.trt = 0, emmeans(nmass, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
nmass.full.70 <- data.frame(n.trt = 70, emmeans(nmass, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
nmass.full.140 <- data.frame(n.trt = 140, emmeans(nmass, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
nmass.full.210 <- data.frame(n.trt = 210, emmeans(nmass, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
nmass.full.280 <- data.frame(n.trt = 280, emmeans(nmass, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
nmass.full.350 <- data.frame(n.trt = 350, emmeans(nmass, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
nmass.full.630 <- data.frame(n.trt = 630, emmeans(nmass, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
nmass.full <- nmass.full.0 %>% full_join(nmass.full.70) %>% 
  full_join(nmass.full.140) %>% full_join(nmass.full.210) %>% 
  full_join(nmass.full.280) %>% full_join(nmass.full.350) %>% 
  full_join(nmass.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

nmass.ntrt.0 <- data.frame(n.trt = 0, emmeans(nmass, ~1, 
                                              at = list(n.trt = 0), type = "response"))
nmass.ntrt.70 <- data.frame(n.trt = 70, emmeans(nmass, ~1, 
                                                at = list(n.trt = 70), type = "response"))
nmass.ntrt.140 <- data.frame(n.trt = 140, emmeans(nmass, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
nmass.ntrt.210 <- data.frame(n.trt = 210, emmeans(nmass, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
nmass.ntrt.280 <- data.frame(n.trt = 280, emmeans(nmass, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
nmass.ntrt.350 <- data.frame(n.trt = 350, emmeans(nmass, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
nmass.ntrt.630 <- data.frame(n.trt = 630, emmeans(nmass, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
nmass.regline <- nmass.ntrt.0 %>% full_join(nmass.ntrt.70) %>%
  full_join(nmass.ntrt.140) %>% full_join(nmass.ntrt.210) %>%
  full_join(nmass.ntrt.280) %>% full_join(nmass.ntrt.350) %>%
  full_join(nmass.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(nmass.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Nmass plot
##########################################################################
nmass.plot <- ggplot(data = df, aes(x = n.trt, y = nmass.focal)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(nmass.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(nmass.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Copy removed outliers and lmer fxn
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(marea)
shapiro.test(residuals(marea))
outlierTest(marea)

test(emtrends(marea, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
marea.full.0 <- data.frame(n.trt = 0, emmeans(marea, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
marea.full.70 <- data.frame(n.trt = 70, emmeans(marea, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
marea.full.140 <- data.frame(n.trt = 140, emmeans(marea, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
marea.full.210 <- data.frame(n.trt = 210, emmeans(marea, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
marea.full.280 <- data.frame(n.trt = 280, emmeans(marea, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
marea.full.350 <- data.frame(n.trt = 350, emmeans(marea, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
marea.full.630 <- data.frame(n.trt = 630, emmeans(marea, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
marea.full <- marea.full.0 %>% full_join(marea.full.70) %>% 
  full_join(marea.full.140) %>% full_join(marea.full.210) %>% 
  full_join(marea.full.280) %>% full_join(marea.full.350) %>% 
  full_join(marea.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

marea.ntrt.0 <- data.frame(n.trt = 0, emmeans(marea, ~1, 
                                              at = list(n.trt = 0), type = "response"))
marea.ntrt.70 <- data.frame(n.trt = 70, emmeans(marea, ~1, 
                                                at = list(n.trt = 70), type = "response"))
marea.ntrt.140 <- data.frame(n.trt = 140, emmeans(marea, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
marea.ntrt.210 <- data.frame(n.trt = 210, emmeans(marea, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
marea.ntrt.280 <- data.frame(n.trt = 280, emmeans(marea, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
marea.ntrt.350 <- data.frame(n.trt = 350, emmeans(marea, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
marea.ntrt.630 <- data.frame(n.trt = 630, emmeans(marea, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
marea.regline <- marea.ntrt.0 %>% full_join(marea.ntrt.70) %>%
  full_join(marea.ntrt.140) %>% full_join(marea.ntrt.210) %>%
  full_join(marea.ntrt.280) %>% full_join(marea.ntrt.350) %>%
  full_join(marea.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(marea.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 30)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linesize = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
marea.plot

##########################################################################
## Chlarea regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$chl.mmolm2[c(111, 113, 114)] <- NA

chlarea <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(chlarea)
shapiro.test(residuals(chlarea))
outlierTest(chlarea)

test(emtrends(chlarea, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
chlarea.full.0 <- data.frame(n.trt = 0, emmeans(chlarea, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
chlarea.full.70 <- data.frame(n.trt = 70, emmeans(chlarea, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
chlarea.full.140 <- data.frame(n.trt = 140, emmeans(chlarea, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
chlarea.full.210 <- data.frame(n.trt = 210, emmeans(chlarea, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
chlarea.full.280 <- data.frame(n.trt = 280, emmeans(chlarea, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
chlarea.full.350 <- data.frame(n.trt = 350, emmeans(chlarea, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
chlarea.full.630 <- data.frame(n.trt = 630, emmeans(chlarea, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
chlarea.full <- chlarea.full.0 %>% full_join(chlarea.full.70) %>% 
  full_join(chlarea.full.140) %>% full_join(chlarea.full.210) %>% 
  full_join(chlarea.full.280) %>% full_join(chlarea.full.350) %>% 
  full_join(chlarea.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

chlarea.ntrt.0 <- data.frame(n.trt = 0, emmeans(chlarea, ~1, 
                                              at = list(n.trt = 0), type = "response"))
chlarea.ntrt.70 <- data.frame(n.trt = 70, emmeans(chlarea, ~1, 
                                                at = list(n.trt = 70), type = "response"))
chlarea.ntrt.140 <- data.frame(n.trt = 140, emmeans(chlarea, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
chlarea.ntrt.210 <- data.frame(n.trt = 210, emmeans(chlarea, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
chlarea.ntrt.280 <- data.frame(n.trt = 280, emmeans(chlarea, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
chlarea.ntrt.350 <- data.frame(n.trt = 350, emmeans(chlarea, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
chlarea.ntrt.630 <- data.frame(n.trt = 630, emmeans(chlarea, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
chlarea.regline <- chlarea.ntrt.0 %>% full_join(chlarea.ntrt.70) %>%
  full_join(chlarea.ntrt.140) %>% full_join(chlarea.ntrt.210) %>%
  full_join(chlarea.ntrt.280) %>% full_join(chlarea.ntrt.350) %>%
  full_join(chlarea.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(chlarea.full) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_inoc" | co2.inoc == "elv_inoc", 
                           "dashed", "solid"))

##########################################################################
## Chl area plot
##########################################################################
chl.plot <- ggplot(data = df, aes(x = n.trt, y = chl.mmolm2)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(chlarea.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(chlarea.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("Chl")["area"]*" (mmol m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linesize = "none") +
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
Anova(vcmax25)
shapiro.test(residuals(vcmax25))
outlierTest(vcmax25)

test(emtrends(vcmax25, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
vcmax.full.0 <- data.frame(n.trt = 0, emmeans(vcmax25, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
vcmax.full.70 <- data.frame(n.trt = 70, emmeans(vcmax25, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
vcmax.full.140 <- data.frame(n.trt = 140, emmeans(vcmax25, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
vcmax.full.210 <- data.frame(n.trt = 210, emmeans(vcmax25, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
vcmax.full.280 <- data.frame(n.trt = 280, emmeans(vcmax25, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
vcmax.full.350 <- data.frame(n.trt = 350, emmeans(vcmax25, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
vcmax.full.630 <- data.frame(n.trt = 630, emmeans(vcmax25, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
vcmax.full <- vcmax.full.0 %>% full_join(vcmax.full.70) %>% 
  full_join(vcmax.full.140) %>% full_join(vcmax.full.210) %>% 
  full_join(vcmax.full.280) %>% full_join(vcmax.full.350) %>% 
  full_join(vcmax.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.ntrt.0 <- data.frame(n.trt = 0, emmeans(vcmax25, ~1, 
                                              at = list(n.trt = 0), type = "response"))
vcmax.ntrt.70 <- data.frame(n.trt = 70, emmeans(vcmax25, ~1, 
                                                at = list(n.trt = 70), type = "response"))
vcmax.ntrt.140 <- data.frame(n.trt = 140, emmeans(vcmax25, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
vcmax.ntrt.210 <- data.frame(n.trt = 210, emmeans(vcmax25, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
vcmax.ntrt.280 <- data.frame(n.trt = 280, emmeans(vcmax25, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
vcmax.ntrt.350 <- data.frame(n.trt = 350, emmeans(vcmax25, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
vcmax.ntrt.630 <- data.frame(n.trt = 630, emmeans(vcmax25, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
vcmax.regline <- vcmax.ntrt.0 %>% full_join(vcmax.ntrt.70) %>%
  full_join(vcmax.ntrt.140) %>% full_join(vcmax.ntrt.210) %>%
  full_join(vcmax.ntrt.280) %>% full_join(vcmax.ntrt.350) %>%
  full_join(vcmax.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
Anova(jmax25)
shapiro.test(residuals(jmax25))
outlierTest(jmax25)

test(emtrends(jmax25, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
jmax.full.0 <- data.frame(n.trt = 0, emmeans(jmax25, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
jmax.full.70 <- data.frame(n.trt = 70, emmeans(jmax25, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
jmax.full.140 <- data.frame(n.trt = 140, emmeans(jmax25, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
jmax.full.210 <- data.frame(n.trt = 210, emmeans(jmax25, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
jmax.full.280 <- data.frame(n.trt = 280, emmeans(jmax25, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
jmax.full.350 <- data.frame(n.trt = 350, emmeans(jmax25, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
jmax.full.630 <- data.frame(n.trt = 630, emmeans(jmax25, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
jmax.full <- jmax.full.0 %>% full_join(jmax.full.70) %>% 
  full_join(jmax.full.140) %>% full_join(jmax.full.210) %>% 
  full_join(jmax.full.280) %>% full_join(jmax.full.350) %>% 
  full_join(jmax.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

jmax.ntrt.0 <- data.frame(n.trt = 0, emmeans(jmax25, ~1, 
                                              at = list(n.trt = 0), type = "response"))
jmax.ntrt.70 <- data.frame(n.trt = 70, emmeans(jmax25, ~1, 
                                                at = list(n.trt = 70), type = "response"))
jmax.ntrt.140 <- data.frame(n.trt = 140, emmeans(jmax25, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
jmax.ntrt.210 <- data.frame(n.trt = 210, emmeans(jmax25, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
jmax.ntrt.280 <- data.frame(n.trt = 280, emmeans(jmax25, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
jmax.ntrt.350 <- data.frame(n.trt = 350, emmeans(jmax25, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
jmax.ntrt.630 <- data.frame(n.trt = 630, emmeans(jmax25, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
jmax.regline <- jmax.ntrt.0 %>% full_join(jmax.ntrt.70) %>%
  full_join(jmax.ntrt.140) %>% full_join(jmax.ntrt.210) %>%
  full_join(jmax.ntrt.280) %>% full_join(jmax.ntrt.350) %>%
  full_join(jmax.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(jmax.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Copy removed outliers and lmer fxn
jvmax <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(jvmax)
shapiro.test(residuals(jvmax))
outlierTest(jmax25)

test(emtrends(jmax25, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
jvmax.full.0 <- data.frame(n.trt = 0, emmeans(jvmax, ~inoc*co2, 
                                             at = list(n.trt = 0), 
                                             type = "response"))
jvmax.full.70 <- data.frame(n.trt = 70, emmeans(jvmax, ~inoc*co2, 
                                               at = list(n.trt = 70), 
                                               type = "response"))
jvmax.full.140 <- data.frame(n.trt = 140, emmeans(jvmax, ~inoc*co2, 
                                                 at = list(n.trt = 140), 
                                                 type = "response")) 
jvmax.full.210 <- data.frame(n.trt = 210, emmeans(jvmax, ~inoc*co2, 
                                                 at = list(n.trt = 210), 
                                                 type = "response"))
jvmax.full.280 <- data.frame(n.trt = 280, emmeans(jvmax, ~inoc*co2, 
                                                 at = list(n.trt = 280), 
                                                 type = "response")) 
jvmax.full.350 <- data.frame(n.trt = 350, emmeans(jvmax, ~inoc*co2, 
                                                 at = list(n.trt = 350), 
                                                 type = "response")) 
jvmax.full.630 <- data.frame(n.trt = 630, emmeans(jvmax, ~inoc*co2, 
                                                 at = list(n.trt = 630),
                                                 type = "response")) 
jvmax.full <- jvmax.full.0 %>% full_join(jvmax.full.70) %>% 
  full_join(jvmax.full.140) %>% full_join(jvmax.full.210) %>% 
  full_join(jvmax.full.280) %>% full_join(jvmax.full.350) %>% 
  full_join(jvmax.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

jvmax.ntrt.0 <- data.frame(n.trt = 0, emmeans(jvmax, ~1, 
                                             at = list(n.trt = 0), type = "response"))
jvmax.ntrt.70 <- data.frame(n.trt = 70, emmeans(jvmax, ~1, 
                                               at = list(n.trt = 70), type = "response"))
jvmax.ntrt.140 <- data.frame(n.trt = 140, emmeans(jvmax, ~1, 
                                                 at = list(n.trt = 140), type = "response"))
jvmax.ntrt.210 <- data.frame(n.trt = 210, emmeans(jvmax, ~1, 
                                                 at = list(n.trt = 210), type = "response"))
jvmax.ntrt.280 <- data.frame(n.trt = 280, emmeans(jvmax, ~1, 
                                                 at = list(n.trt = 280), type = "response"))
jvmax.ntrt.350 <- data.frame(n.trt = 350, emmeans(jvmax, ~1, 
                                                 at = list(n.trt = 350), type = "response"))
jvmax.ntrt.630 <- data.frame(n.trt = 630, emmeans(jvmax, ~1, 
                                                 at = list(n.trt = 630), type = "response"))
jvmax.regline <- jmax.ntrt.0 %>% full_join(jvmax.ntrt.70) %>%
  full_join(jvmax.ntrt.140) %>% full_join(jvmax.ntrt.210) %>%
  full_join(jvmax.ntrt.280) %>% full_join(jvmax.ntrt.350) %>%
  full_join(jvmax.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(jvmax.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Prop leaf N in photosynthesis  regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$p.photo[50] <- NA

p.photo <- lmer(p.photo ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(p.photo)
shapiro.test(residuals(p.photo))
outlierTest(p.photo)

test(emtrends(p.photo, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
p.photo.full.0 <- data.frame(n.trt = 0, emmeans(p.photo, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
p.photo.full.70 <- data.frame(n.trt = 70, emmeans(p.photo, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
p.photo.full.140 <- data.frame(n.trt = 140, emmeans(p.photo, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
p.photo.full.210 <- data.frame(n.trt = 210, emmeans(p.photo, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
p.photo.full.280 <- data.frame(n.trt = 280, emmeans(p.photo, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
p.photo.full.350 <- data.frame(n.trt = 350, emmeans(p.photo, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
p.photo.full.630 <- data.frame(n.trt = 630, emmeans(p.photo, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
p.photo.full <- p.photo.full.0 %>% full_join(p.photo.full.70) %>% 
  full_join(p.photo.full.140) %>% full_join(p.photo.full.210) %>% 
  full_join(p.photo.full.280) %>% full_join(p.photo.full.350) %>% 
  full_join(p.photo.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.photo.ntrt.0 <- data.frame(n.trt = 0, emmeans(p.photo, ~1, 
                                              at = list(n.trt = 0), type = "response"))
p.photo.ntrt.70 <- data.frame(n.trt = 70, emmeans(p.photo, ~1, 
                                                at = list(n.trt = 70), type = "response"))
p.photo.ntrt.140 <- data.frame(n.trt = 140, emmeans(p.photo, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
p.photo.ntrt.210 <- data.frame(n.trt = 210, emmeans(p.photo, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
p.photo.ntrt.280 <- data.frame(n.trt = 280, emmeans(p.photo, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
p.photo.ntrt.350 <- data.frame(n.trt = 350, emmeans(p.photo, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
p.photo.ntrt.630 <- data.frame(n.trt = 630, emmeans(p.photo, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
p.photo.regline <- p.photo.ntrt.0 %>% full_join(p.photo.ntrt.70) %>%
  full_join(p.photo.ntrt.140) %>% full_join(p.photo.ntrt.210) %>%
  full_join(p.photo.ntrt.280) %>% full_join(p.photo.ntrt.350) %>%
  full_join(p.photo.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.photo.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_no.inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Copy removed outliers and lmer fxn
df$p.rubisco[50] <- NA

p.rub <- lmer(p.rubisco ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(p.rub)
shapiro.test(residuals(p.rub))
outlierTest(p.rub)

test(emtrends(p.rub, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
p.rub.full.0 <- data.frame(n.trt = 0, emmeans(p.rub, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
p.rub.full.70 <- data.frame(n.trt = 70, emmeans(p.rub, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
p.rub.full.140 <- data.frame(n.trt = 140, emmeans(p.rub, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
p.rub.full.210 <- data.frame(n.trt = 210, emmeans(p.rub, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
p.rub.full.280 <- data.frame(n.trt = 280, emmeans(p.rub, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
p.rub.full.350 <- data.frame(n.trt = 350, emmeans(p.rub, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
p.rub.full.630 <- data.frame(n.trt = 630, emmeans(p.rub, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
p.rub.full <- p.rub.full.0 %>% full_join(p.rub.full.70) %>% 
  full_join(p.rub.full.140) %>% full_join(p.rub.full.210) %>% 
  full_join(p.rub.full.280) %>% full_join(p.rub.full.350) %>% 
  full_join(p.rub.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.rub.ntrt.0 <- data.frame(n.trt = 0, emmeans(p.rub, ~1, 
                                              at = list(n.trt = 0), type = "response"))
p.rub.ntrt.70 <- data.frame(n.trt = 70, emmeans(p.rub, ~1, 
                                                at = list(n.trt = 70), type = "response"))
p.rub.ntrt.140 <- data.frame(n.trt = 140, emmeans(p.rub, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
p.rub.ntrt.210 <- data.frame(n.trt = 210, emmeans(p.rub, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
p.rub.ntrt.280 <- data.frame(n.trt = 280, emmeans(p.rub, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
p.rub.ntrt.350 <- data.frame(n.trt = 350, emmeans(p.rub, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
p.rub.ntrt.630 <- data.frame(n.trt = 630, emmeans(p.rub, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
p.rub.regline <- p.rub.ntrt.0 %>% full_join(p.rub.ntrt.70) %>%
  full_join(p.rub.ntrt.140) %>% full_join(p.rub.ntrt.210) %>%
  full_join(p.rub.ntrt.280) %>% full_join(p.rub.ntrt.350) %>%
  full_join(p.rub.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.rub.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_no.inoc" | co2.inoc == "amb_no.inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Prop leaf N in structure regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$p.structure[50] <- NA

p.str <- lmer(p.structure ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(p.str)
shapiro.test(residuals(p.str))
outlierTest(p.str)

test(emtrends(p.str, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
p.str.full.0 <- data.frame(n.trt = 0, emmeans(p.str, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
p.str.full.70 <- data.frame(n.trt = 70, emmeans(p.str, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
p.str.full.140 <- data.frame(n.trt = 140, emmeans(p.str, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
p.str.full.210 <- data.frame(n.trt = 210, emmeans(p.str, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
p.str.full.280 <- data.frame(n.trt = 280, emmeans(p.str, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
p.str.full.350 <- data.frame(n.trt = 350, emmeans(p.str, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
p.str.full.630 <- data.frame(n.trt = 630, emmeans(p.str, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
p.str.full <- p.str.full.0 %>% full_join(p.str.full.70) %>% 
  full_join(p.str.full.140) %>% full_join(p.str.full.210) %>% 
  full_join(p.str.full.280) %>% full_join(p.str.full.350) %>% 
  full_join(p.str.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

p.str.ntrt.0 <- data.frame(n.trt = 0, emmeans(p.str, ~1, 
                                              at = list(n.trt = 0), type = "response"))
p.str.ntrt.70 <- data.frame(n.trt = 70, emmeans(p.str, ~1, 
                                                at = list(n.trt = 70), type = "response"))
p.str.ntrt.140 <- data.frame(n.trt = 140, emmeans(p.str, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
p.str.ntrt.210 <- data.frame(n.trt = 210, emmeans(p.str, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
p.str.ntrt.280 <- data.frame(n.trt = 280, emmeans(p.str, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
p.str.ntrt.350 <- data.frame(n.trt = 350, emmeans(p.str, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
p.str.ntrt.630 <- data.frame(n.trt = 630, emmeans(p.str, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
p.str.regline <- p.str.ntrt.0 %>% full_join(p.str.ntrt.70) %>%
  full_join(p.str.ntrt.140) %>% full_join(p.str.ntrt.210) %>%
  full_join(p.str.ntrt.280) %>% full_join(p.str.ntrt.350) %>%
  full_join(p.str.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(p.str.full) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Prop leaf N in structure plot
##########################################################################
p.str.plot <- ggplot(data = df, aes(x = n.trt, y = p.structure)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(p.str.regline, co2.inoc != "overall"),
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(p.str.regline, co2.inoc != "overall"),
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.06)) +
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
df$pnue[50] <- NA

pnue <- lmer(pnue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(pnue)
shapiro.test(residuals(pnue))
outlierTest(pnue)

test(emtrends(pnue, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
pnue.full.0 <- data.frame(n.trt = 0, emmeans(pnue, ~inoc*co2, 
                                              at = list(n.trt = 0), 
                                              type = "response"))
pnue.full.70 <- data.frame(n.trt = 70, emmeans(pnue, ~inoc*co2, 
                                                at = list(n.trt = 70), 
                                                type = "response"))
pnue.full.140 <- data.frame(n.trt = 140, emmeans(pnue, ~inoc*co2, 
                                                  at = list(n.trt = 140), 
                                                  type = "response")) 
pnue.full.210 <- data.frame(n.trt = 210, emmeans(pnue, ~inoc*co2, 
                                                  at = list(n.trt = 210), 
                                                  type = "response"))
pnue.full.280 <- data.frame(n.trt = 280, emmeans(pnue, ~inoc*co2, 
                                                  at = list(n.trt = 280), 
                                                  type = "response")) 
pnue.full.350 <- data.frame(n.trt = 350, emmeans(pnue, ~inoc*co2, 
                                                  at = list(n.trt = 350), 
                                                  type = "response")) 
pnue.full.630 <- data.frame(n.trt = 630, emmeans(pnue, ~inoc*co2, 
                                                  at = list(n.trt = 630),
                                                  type = "response")) 
pnue.full <- pnue.full.0 %>% full_join(pnue.full.70) %>% 
  full_join(pnue.full.140) %>% full_join(pnue.full.210) %>% 
  full_join(pnue.full.280) %>% full_join(pnue.full.350) %>% 
  full_join(pnue.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

pnue.ntrt.0 <- data.frame(n.trt = 0, emmeans(pnue, ~1, 
                                              at = list(n.trt = 0), type = "response"))
pnue.ntrt.70 <- data.frame(n.trt = 70, emmeans(pnue, ~1, 
                                                at = list(n.trt = 70), type = "response"))
pnue.ntrt.140 <- data.frame(n.trt = 140, emmeans(pnue, ~1, 
                                                  at = list(n.trt = 140), type = "response"))
pnue.ntrt.210 <- data.frame(n.trt = 210, emmeans(pnue, ~1, 
                                                  at = list(n.trt = 210), type = "response"))
pnue.ntrt.280 <- data.frame(n.trt = 280, emmeans(pnue, ~1, 
                                                  at = list(n.trt = 280), type = "response"))
pnue.ntrt.350 <- data.frame(n.trt = 350, emmeans(pnue, ~1, 
                                                  at = list(n.trt = 350), type = "response"))
pnue.ntrt.630 <- data.frame(n.trt = 630, emmeans(pnue, ~1, 
                                                  at = list(n.trt = 630), type = "response"))
pnue.regline <- pnue.ntrt.0 %>% full_join(pnue.ntrt.70) %>%
  full_join(pnue.ntrt.140) %>% full_join(pnue.ntrt.210) %>%
  full_join(pnue.ntrt.280) %>% full_join(pnue.ntrt.350) %>%
  full_join(pnue.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(pnue.full) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc" | co2.inoc == "elv_no.inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
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
## Copy removed outliers and lmer fxn
df$pnue[50] <- NA

iwue <- lmer(iwue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(iwue)
shapiro.test(residuals(iwue))
outlierTest(iwue)

test(emtrends(iwue, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
iwue.full.0 <- data.frame(n.trt = 0, emmeans(iwue, ~inoc*co2, 
                                             at = list(n.trt = 0), 
                                             type = "response"))
iwue.full.70 <- data.frame(n.trt = 70, emmeans(iwue, ~inoc*co2, 
                                               at = list(n.trt = 70), 
                                               type = "response"))
iwue.full.140 <- data.frame(n.trt = 140, emmeans(iwue, ~inoc*co2, 
                                                 at = list(n.trt = 140), 
                                                 type = "response")) 
iwue.full.210 <- data.frame(n.trt = 210, emmeans(iwue, ~inoc*co2, 
                                                 at = list(n.trt = 210), 
                                                 type = "response"))
iwue.full.280 <- data.frame(n.trt = 280, emmeans(iwue, ~inoc*co2, 
                                                 at = list(n.trt = 280), 
                                                 type = "response")) 
iwue.full.350 <- data.frame(n.trt = 350, emmeans(iwue, ~inoc*co2, 
                                                 at = list(n.trt = 350), 
                                                 type = "response")) 
iwue.full.630 <- data.frame(n.trt = 630, emmeans(iwue, ~inoc*co2, 
                                                 at = list(n.trt = 630),
                                                 type = "response")) 
iwue.full <- iwue.full.0 %>% full_join(iwue.full.70) %>% 
  full_join(iwue.full.140) %>% full_join(iwue.full.210) %>% 
  full_join(iwue.full.280) %>% full_join(iwue.full.350) %>% 
  full_join(iwue.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

iwue.ntrt.0 <- data.frame(n.trt = 0, emmeans(iwue, ~1, 
                                             at = list(n.trt = 0), type = "response"))
iwue.ntrt.70 <- data.frame(n.trt = 70, emmeans(iwue, ~1, 
                                               at = list(n.trt = 70), type = "response"))
iwue.ntrt.140 <- data.frame(n.trt = 140, emmeans(iwue, ~1, 
                                                 at = list(n.trt = 140), type = "response"))
iwue.ntrt.210 <- data.frame(n.trt = 210, emmeans(iwue, ~1, 
                                                 at = list(n.trt = 210), type = "response"))
iwue.ntrt.280 <- data.frame(n.trt = 280, emmeans(iwue, ~1, 
                                                 at = list(n.trt = 280), type = "response"))
iwue.ntrt.350 <- data.frame(n.trt = 350, emmeans(iwue, ~1, 
                                                 at = list(n.trt = 350), type = "response"))
iwue.ntrt.630 <- data.frame(n.trt = 630, emmeans(iwue, ~1, 
                                                 at = list(n.trt = 630), type = "response"))
iwue.regline <- iwue.ntrt.0 %>% full_join(iwue.ntrt.70) %>%
  full_join(iwue.ntrt.140) %>% full_join(iwue.ntrt.210) %>%
  full_join(iwue.ntrt.280) %>% full_join(iwue.ntrt.350) %>%
  full_join(iwue.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(iwue.full) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc" | co2.inoc == "elv_no.inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
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
        panel.border = element_rect(size = 1.25))
iwue.plot

##########################################################################
## Narea:gs regression line prep
##########################################################################
narea.gs <- lmer(log(narea.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
Anova(narea.gs)
shapiro.test(residuals(narea.gs))
outlierTest(narea.gs)

test(emtrends(narea.gs, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
narea.gs.full.0 <- data.frame(n.trt = 0, emmeans(narea.gs, ~inoc*co2, 
                                             at = list(n.trt = 0), 
                                             type = "response"))
narea.gs.full.70 <- data.frame(n.trt = 70, emmeans(narea.gs, ~inoc*co2, 
                                               at = list(n.trt = 70), 
                                               type = "response"))
narea.gs.full.140 <- data.frame(n.trt = 140, emmeans(narea.gs, ~inoc*co2, 
                                                 at = list(n.trt = 140), 
                                                 type = "response")) 
narea.gs.full.210 <- data.frame(n.trt = 210, emmeans(narea.gs, ~inoc*co2, 
                                                 at = list(n.trt = 210), 
                                                 type = "response"))
narea.gs.full.280 <- data.frame(n.trt = 280, emmeans(narea.gs, ~inoc*co2, 
                                                 at = list(n.trt = 280), 
                                                 type = "response")) 
narea.gs.full.350 <- data.frame(n.trt = 350, emmeans(narea.gs, ~inoc*co2, 
                                                 at = list(n.trt = 350), 
                                                 type = "response")) 
narea.gs.full.630 <- data.frame(n.trt = 630, emmeans(narea.gs, ~inoc*co2, 
                                                 at = list(n.trt = 630),
                                                 type = "response")) 
narea.gs.full <- narea.gs.full.0 %>% full_join(narea.gs.full.70) %>% 
  full_join(narea.gs.full.140) %>% full_join(narea.gs.full.210) %>% 
  full_join(narea.gs.full.280) %>% full_join(narea.gs.full.350) %>% 
  full_join(narea.gs.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.gs.ntrt.0 <- data.frame(n.trt = 0, emmeans(narea.gs, ~1, 
                                             at = list(n.trt = 0), type = "response"))
narea.gs.ntrt.70 <- data.frame(n.trt = 70, emmeans(narea.gs, ~1, 
                                               at = list(n.trt = 70), type = "response"))
narea.gs.ntrt.140 <- data.frame(n.trt = 140, emmeans(narea.gs, ~1, 
                                                 at = list(n.trt = 140), type = "response"))
narea.gs.ntrt.210 <- data.frame(n.trt = 210, emmeans(narea.gs, ~1, 
                                                 at = list(n.trt = 210), type = "response"))
narea.gs.ntrt.280 <- data.frame(n.trt = 280, emmeans(narea.gs, ~1, 
                                                 at = list(n.trt = 280), type = "response"))
narea.gs.ntrt.350 <- data.frame(n.trt = 350, emmeans(narea.gs, ~1, 
                                                 at = list(n.trt = 350), type = "response"))
narea.gs.ntrt.630 <- data.frame(n.trt = 630, emmeans(narea.gs, ~1, 
                                                 at = list(n.trt = 630), type = "response"))
narea.gs.regline <- narea.gs.ntrt.0 %>% full_join(narea.gs.ntrt.70) %>%
  full_join(narea.gs.ntrt.140) %>% full_join(narea.gs.ntrt.210) %>%
  full_join(narea.gs.ntrt.280) %>% full_join(narea.gs.ntrt.350) %>%
  full_join(narea.gs.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.gs.full) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc" | co2.inoc == "elv_no.inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = cbbPalette3,
                     labels = c("Ambient, inoculated",
                                "Ambient, not inoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 32, 8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*": "*italic("g")["s"]*" ()")),
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
Anova(vcmax.gs)
shapiro.test(residuals(vcmax.gs))
outlierTest(vcmax.gs)

test(emtrends(vcmax.gs, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax.gs.full.0 <- data.frame(n.trt = 0, emmeans(vcmax.gs, ~inoc*co2, 
                                                 at = list(n.trt = 0), 
                                                 type = "response"))
vcmax.gs.full.70 <- data.frame(n.trt = 70, emmeans(vcmax.gs, ~inoc*co2, 
                                                   at = list(n.trt = 70), 
                                                   type = "response"))
vcmax.gs.full.140 <- data.frame(n.trt = 140, emmeans(vcmax.gs, ~inoc*co2, 
                                                     at = list(n.trt = 140), 
                                                     type = "response")) 
vcmax.gs.full.210 <- data.frame(n.trt = 210, emmeans(vcmax.gs, ~inoc*co2, 
                                                     at = list(n.trt = 210), 
                                                     type = "response"))
vcmax.gs.full.280 <- data.frame(n.trt = 280, emmeans(vcmax.gs, ~inoc*co2, 
                                                     at = list(n.trt = 280), 
                                                     type = "response")) 
vcmax.gs.full.350 <- data.frame(n.trt = 350, emmeans(vcmax.gs, ~inoc*co2, 
                                                     at = list(n.trt = 350), 
                                                     type = "response")) 
vcmax.gs.full.630 <- data.frame(n.trt = 630, emmeans(vcmax.gs, ~inoc*co2, 
                                                     at = list(n.trt = 630),
                                                     type = "response")) 
vcmax.gs.full <- vcmax.gs.full.0 %>% full_join(vcmax.gs.full.70) %>% 
  full_join(vcmax.gs.full.140) %>% full_join(vcmax.gs.full.210) %>% 
  full_join(vcmax.gs.full.280) %>% full_join(vcmax.gs.full.350) %>% 
  full_join(vcmax.gs.full.630) %>% 
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.gs.ntrt.0 <- data.frame(n.trt = 0, emmeans(vcmax.gs, ~1, 
                                                 at = list(n.trt = 0), type = "response"))
vcmax.gs.ntrt.70 <- data.frame(n.trt = 70, emmeans(vcmax.gs, ~1, 
                                                   at = list(n.trt = 70), type = "response"))
vcmax.gs.ntrt.140 <- data.frame(n.trt = 140, emmeans(vcmax.gs, ~1, 
                                                     at = list(n.trt = 140), type = "response"))
vcmax.gs.ntrt.210 <- data.frame(n.trt = 210, emmeans(vcmax.gs, ~1, 
                                                     at = list(n.trt = 210), type = "response"))
vcmax.gs.ntrt.280 <- data.frame(n.trt = 280, emmeans(vcmax.gs, ~1, 
                                                     at = list(n.trt = 280), type = "response"))
vcmax.gs.ntrt.350 <- data.frame(n.trt = 350, emmeans(vcmax.gs, ~1, 
                                                     at = list(n.trt = 350), type = "response"))
vcmax.gs.ntrt.630 <- data.frame(n.trt = 630, emmeans(vcmax.gs, ~1, 
                                                     at = list(n.trt = 630), type = "response"))
vcmax.gs.regline <- vcmax.gs.ntrt.0 %>% full_join(vcmax.gs.ntrt.70) %>%
  full_join(vcmax.gs.ntrt.140) %>% full_join(vcmax.gs.ntrt.210) %>%
  full_join(vcmax.gs.ntrt.280) %>% full_join(vcmax.gs.ntrt.350) %>%
  full_join(vcmax.gs.ntrt.630) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.gs.full) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc" | co2.inoc == "elv_no.inoc", 
                           "dashed", "solid")) %>%
  select(n.trt, co2, inoc, co2.inoc, everything(), -X1)

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
                               "Ambient, not inoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 200)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*": "*italic("g")["s"]*" (μmol CO"["2"]*" mol"^"-1"*"H"["2"]*"O)")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
vcmax.gs.plot

##########################################################################
## Figure 1: whole plant plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig1_wholePlant.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(ncost.plot, tla.plot, tbio.plot, ncol = 2, nrow = 2,
          common.legend = TRUE, align = "hv",
          legend = "right")
dev.off()

##########################################################################
## Figure 2: leaf N plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig2_leafN.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(narea.plot, nmass.plot, marea.plot, chl.plot,
          ncol = 2, nrow = 2,
          common.legend = TRUE, align = "hv",
          legend = "right")
dev.off()

##########################################################################
## Figure 3: leaf physiology plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig3_photo.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(vcmax.plot, jmax.plot, jvmax.plot, ncol = 2, nrow = 2,
          common.legend = TRUE, align = "hv",
          legend = "right")
dev.off()

##########################################################################
## Figure 4: propN photosynthesis/structure
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig4_propN.png",
    height = 4, width = 12, units = "in", res = 600)
ggarrange(p.photo.plot, p.str.plot, ncol = 2, nrow = 1,
          common.legend = TRUE, align = "hv",
          legend = "right")
dev.off()

##########################################################################
## Figure 5: PNUE/iWUE tradeoffs
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig5_PNUE_iWUE.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(pnue.plot, iwue.plot, narea.gs.plot, vcmax.gs.plot,
          ncol = 2, nrow = 2, common.legend = TRUE, align = "hv",
          legend = "right")
dev.off()




