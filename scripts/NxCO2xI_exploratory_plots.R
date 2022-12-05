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


df$ncost[df$ncost > 20] <- NA

##########################################################################
## Ncost:Vcmax25 regline prep
##########################################################################
vcmax.ncost <- lmer(vcmax25 ~ ncost * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(vcmax.ncost))
test(emtrends(vcmax.ncost, ~co2*inoc, "ncost"))

## Emmean fxns for regression lines + error ribbons
vcmax.ncost.full <- data.frame(emmeans(vcmax.ncost, ~inoc*co2, "ncost",
                                       at = list(ncost = seq(0, 20, 0.1)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.ncost.regline <- data.frame(emmeans(vcmax.ncost, ~1, "ncost",
                                          at = list(ncost = seq(0, 20, 0.1)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.ncost.full) %>%
  dplyr::select(ncost, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Ncost:Vcmax25 plot
##########################################################################
vcmax.ncost.plot <- ggplot(data = df, 
                           aes(x = ncost, y = vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(vcmax.ncost.regline, co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(vcmax.ncost.regline, co2.inoc == "overall"),
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
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
  scale_y_continuous(limits = c(-30, 150), breaks = seq(-30, 150, 30)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("V")["cmax25"]*" (μmol m"^"-2"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
vcmax.ncost.plot

##########################################################################
## Rd.vcmax:Vcmax25 regline prep
##########################################################################
vcmax.rd.vcmax <- lmer(vcmax25 ~ rd25.vcmax25 * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(vcmax.rd.vcmax))
test(emtrends(vcmax.rd.vcmax, ~co2*inoc, "rd25.vcmax25"))

## Emmean fxns for regression lines + error ribbons
vcmax.rd.vcmax.full <- data.frame(emmeans(vcmax.rd.vcmax, ~inoc*co2, "rd25.vcmax25",
                                       at = list(rd25.vcmax25 = seq(0, 0.09, 0.001)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.rd.vcmax.regline <- data.frame(emmeans(vcmax.rd.vcmax, ~1, "rd25.vcmax25",
                                          at = list(rd25.vcmax25 = seq(0, 0.09, 0.001)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.rd.vcmax.full) %>%
  dplyr::select(rd25.vcmax25, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc", 
                           "dashed", "solid"))

##########################################################################
## Ncost:Vcmax25 plot
##########################################################################
vcmax.ncost.plot <- ggplot(data = df, 
                           aes(x = rd25.vcmax25, y = vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(vcmax.rd.vcmax.regline, co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(vcmax.rd.vcmax.regline, co2.inoc == "overall"),
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
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
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 30)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("V")["cmax25"]*" (μmol m"^"-2"*"s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
vcmax.ncost.plot

##########################################################################
## Ncost:Narea regline prep
##########################################################################

narea.ncost <- lmer(sqrt(narea) ~ ncost * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(narea.ncost))
test(emtrends(narea.ncost, ~co2*inoc, "ncost"))

## Emmean fxns for regression lines + error ribbons
narea.ncost.full <- data.frame(emmeans(narea.ncost, ~inoc*co2, "ncost",
                                       at = list(ncost = seq(0, 20, 0.1)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.ncost.regline <- data.frame(emmeans(narea.ncost, ~1, "ncost",
                                          at = list(ncost = seq(0, 20, 0.1)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.ncost.full) %>%
  dplyr::select(ncost, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Ncost:Narea plot
##########################################################################
narea.ncost.plot <- ggplot(data = df, 
                           aes(x = ncost, y = narea)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(narea.ncost.regline, co2.inoc == "overall"),
              aes(y = response), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(narea.ncost.regline, co2.inoc == "overall"),
              aes(y = response, ymin = lower.CL, ymax = upper.CL), 
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
  scale_y_continuous(limits = c(0, 3.2), breaks = seq(0, 3, 1)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
narea.ncost.plot

##########################################################################
## Rd.vcmax:Narea regline prep
##########################################################################
narea.rd.vcmax <- lmer(narea ~ rd25.vcmax25 * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(narea.rd.vcmax))
test(emtrends(narea.rd.vcmax, ~co2*inoc, "rd25.vcmax25"))

## Emmean fxns for regression lines + error ribbons
narea.rd.vcmax.full <- data.frame(emmeans(narea.rd.vcmax, ~inoc*co2, "rd25.vcmax25",
                                       at = list(rd25.vcmax25 = seq(0, 0.09, 0.001)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.rd.vcmax.regline <- data.frame(emmeans(narea.rd.vcmax, ~1, "rd25.vcmax25",
                                          at = list(rd25.vcmax25 = seq(0, 0.09, 0.001)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.rd.vcmax.full) %>%
  dplyr::select(rd25.vcmax25, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Rd.vcmax:Narea plot
##########################################################################
narea.rd.vcmax.plot <- ggplot(data = df, 
                           aes(x = rd25.vcmax25, y = narea)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(narea.rd.vcmax.regline, co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(narea.rd.vcmax.regline, co2.inoc == "overall"),
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
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
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
narea.rd.vcmax.plot

##########################################################################
## Ncost:Narea regline prep
##########################################################################
jmax25.vcmax25.ncost <- lmer(log(jmax25.vcmax25) ~ ncost * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(jmax25.vcmax25.ncost))
test(emtrends(jmax25.vcmax25.ncost, ~co2*inoc, "ncost"))

## Emmean fxns for regression lines + error ribbons
jmax25.vcmax25.ncost.full <- data.frame(emmeans(jmax25.vcmax25.ncost, ~inoc*co2, "ncost",
                                       at = list(ncost = seq(0, 20, 0.1)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

jmax25.vcmax25.ncost.regline <- data.frame(emmeans(jmax25.vcmax25.ncost, ~1, "ncost",
                                          at = list(ncost = seq(0, 20, 0.1)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(jmax25.vcmax25.ncost.full) %>%
  dplyr::select(ncost, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb.inoc",
                           "dashed", "solid"))

##########################################################################
## Ncost:Narea plot
##########################################################################

jmax25.vcmax25.ncost.plot <- ggplot(data = df, 
                           aes(x = ncost, y = jmax25.vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(jmax25.vcmax25.ncost.regline, co2.inoc == "overall"),
              aes(y = response), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(jmax25.vcmax25.ncost.regline, co2.inoc == "overall"),
              aes(y = response, ymin = lower.CL, ymax = upper.CL), 
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
  scale_y_continuous(limits = c(1.5, 2.3), breaks = seq(1.5, 2.25, 0.25)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("J")["max25"]*":"*italic(V)["cmax25"])),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
jmax25.vcmax25.ncost.plot

##########################################################################
## Ncost:Rd25.vcmax25 regline prep
##########################################################################
jmax25.vcmax25.rd25.vcmax25 <- lmer(jmax25.vcmax25 ~ rd25.vcmax25 * 
                                      co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(jmax25.vcmax25.rd25.vcmax25))
test(emtrends(jmax25.vcmax25.rd25.vcmax25, ~co2*inoc, "rd25.vcmax25"))

## Emmean fxns for regression lines + error ribbons
jmax25.vcmax25.rd25.vcmax25.full <- data.frame(
  emmeans(jmax25.vcmax25.rd25.vcmax25, 
          ~inoc*co2, "rd25.vcmax25",
          at = list(rd25.vcmax25 = seq(0, 0.09, 0.001)),
          type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

jmax25.vcmax25.rd25.vcmax25.regline <- data.frame(
  emmeans(jmax25.vcmax25.rd25.vcmax25, ~1, "rd25.vcmax25",
          at = list(rd25.vcmax25 = seq(0, 0.09, 0.001)),
          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(jmax25.vcmax25.rd25.vcmax25.full) %>%
  dplyr::select(rd25.vcmax25, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb.inoc",
                           "dashed", "solid"))

##########################################################################
## Ncost:Narea plot
##########################################################################

jmax25.vcmax25.rd25.vcmax25.plot <- ggplot(data = df, 
                                           aes(x = rd25.vcmax25, 
                                               y = jmax25.vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(jmax25.vcmax25.rd25.vcmax25.regline, 
                            co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(jmax25.vcmax25.rd25.vcmax25.regline, 
                            co2.inoc == "overall"),
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
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
  scale_y_continuous(limits = c(1.5, 2.3), breaks = seq(1.5, 2.25, 0.25)) +
  scale_x_continuous(limits = c(0, 0.09), breaks = seq(0, 0.09, 0.03)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("J")["max25"]*":"*italic(V)["cmax25"])),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
jmax25.vcmax25.rd25.vcmax25.plot

##########################################################################
## Ncost:Rd.vcmax regline prep
##########################################################################
df$rd25.vcmax25 <- df$rd25 / df$vcmax25
df$rd25.vcmax25[df$rd25.vcmax25 < 0] <- NA
df$rd25.vcmax25[c(39, 40, 42)] <- NA

rd25.vcmax25.ncost <- lmer(rd25.vcmax25 ~ ncost * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(rd25.vcmax25.ncost))
test(emtrends(rd25.vcmax25.ncost, ~co2*inoc, "ncost"))

## Emmean fxns for regression lines + error ribbons
rd25.vcmax25.ncost.full <- data.frame(emmeans(rd25.vcmax25.ncost, ~inoc*co2, "ncost",
                                                at = list(ncost = seq(0, 20, 0.1)),
                                                type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

rd25.vcmax25.ncost.regline <- data.frame(emmeans(rd25.vcmax25.ncost, ~1, "ncost",
                                                   at = list(ncost = seq(0, 20, 0.1)),
                                                   type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(rd25.vcmax25.ncost.full) %>%
  dplyr::select(ncost, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb.inoc",
                           "dashed", "solid"))


##########################################################################
## Ncost:Rd.vcmax plot
##########################################################################
rd25.vcmax25.ncost.plot <-ggplot(data = df, aes(x = ncost, y = rd25.vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(rd25.vcmax25.ncost.regline, co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(rd25.vcmax25.ncost.regline, co2.inoc == "overall"),
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = cbbPalette3,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(0, 0.09), breaks = seq(0, 0.09, 0.03)) +
  labs(x = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       y = expression(bold(italic("R")["d25"]*":"*italic(V)["cmax25"])),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
rd25.vcmax25.ncost.plot

##########################################################################
## Exploratory figs 
##########################################################################

png("../working_drafts/figs/NxCO2xI_expl_rdvcmax_ncost.png",
    width = 9, height = 5.5, units = "in", res = 600)
rd25.vcmax25.ncost.plot
dev.off()

png("../working_drafts/figs/NxCO2xI_expl_figs.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(narea.ncost.plot, vcmax.ncost.plot, jmax25.vcmax25.ncost.plot,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "right",
          align = "hv")
dev.off()



