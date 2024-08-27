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
## Ncost:Narea regline prep
##########################################################################
narea.beta <- lmer(narea ~ chi * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(narea.beta))
test(emtrends(narea.beta, ~co2*inoc, "chi"))

## Emmean fxns for regression lines + error ribbons
narea.beta.full <- data.frame(emmeans(narea.beta, ~inoc*co2, "chi",
                                       at = list(chi = seq(0.4, 0.8, 0.01)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

narea.beta.regline <- data.frame(emmeans(narea.beta, ~1, "chi",
                                          at = list(chi = seq(0.4, 0.8, 0.01)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(narea.beta.full) %>%
  dplyr::select(chi, co2, inoc, co2.inoc, everything(), -X1)

##########################################################################
## Ncost:Narea plot
##########################################################################
narea.beta.plot <- ggplot(data = df, 
                           aes(x = beta, y = ncost)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(narea.beta.regline, co2.inoc == "overall"),
              aes(x = chi, y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(narea.beta.regline, co2.inoc == "overall"),
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
## Ncost:Jvmax25 regline prep
##########################################################################
jmax25.vcmax25.ncost <- lmer(sqrt(jmax25.vcmax25) ~ ncost * co2 * inoc + (1|rack:co2), data = df)
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
## Ncost: Jvmax25 plot
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
## Ncost:Narea regline prep
##########################################################################
ncost.beta <- lmer(ncost ~ log(beta) * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(ncost.beta))
test(emtrends(ncost.beta, ~co2*inoc, "beta"))

## Emmean fxns for regression lines + error ribbons
ncost.beta.full <- data.frame(emmeans(ncost.beta, ~inoc*co2, "beta",
                                      at = list(beta = seq(0, 400, 1)),
                                      type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

ncost.beta.elv_inoc <- subset(ncost.beta.full, co2.inoc == "elv_inoc" & beta < 150)
ncost.beta.elv_no.inoc <- subset(ncost.beta.full, co2.inoc == "elv_no.inoc" & beta < 100)
ncost.beta.amb_inoc <- subset(ncost.beta.full, co2.inoc == "amb_inoc" & beta < 225)
ncost.beta.amb_no.inoc <- subset(ncost.beta.full, co2.inoc == "amb_no.inoc" & beta > 100)

ncost.beta.cleaned <- ncost.beta.elv_inoc %>% full_join(ncost.beta.elv_no.inoc) %>%
  full_join(ncost.beta.amb_inoc) %>% full_join(ncost.beta.amb_no.inoc)


ncost.beta.regline <- data.frame(emmeans(ncost.beta, ~1, "beta",
                                         at = list(beta = seq(0, 400, 1)),
                                         type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(ncost.beta.cleaned) %>%
  dplyr::select(beta, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "elv_inoc" | co2.inoc == "amb_inoc",
                           "dashed", "solid"))

##########################################################################
## Ncost:Narea plot
##########################################################################
ncost.beta.plot <- ggplot(data = df, 
                          aes(x = beta, y = ncost)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(ncost.beta.regline, co2.inoc != "overall"),
              aes(x = beta, y = response, color = co2.inoc, linetype = linetype),
              size = 1.5, se = FALSE) +
  geom_ribbon(data = subset(ncost.beta.regline, co2.inoc != "overall"),
              aes(y = response, ymin = lower.CL, ymax = upper.CL, fill = co2.inoc), 
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
  scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 8)) +
  labs(x = expression(beta),
       y = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))

png("../working_drafts/figs/NxCO2xI_beta_ncost.png",
    width = 8, height = 5, units = "in", res = 600)
ncost.beta.plot
dev.off()

##########################################################################
## beta:Vcmax25 regline prep
##########################################################################
vcmax.beta <- lmer(vcmax25 ~ beta * co2 * inoc + (1|rack:co2), data = df)
shapiro.test(residuals(vcmax.beta))
test(emtrends(vcmax.beta, ~co2*inoc, "beta"))

## Emmean fxns for regression lines + error ribbons
vcmax.beta.full <- data.frame(emmeans(vcmax.beta, ~inoc*co2, "beta",
                                       at = list(beta = seq(0, 400, 1)),
                                       type = "response")) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE)

vcmax.beta.regline <- data.frame(emmeans(vcmax.beta, ~1, "beta",
                                          at = list(beta = seq(0, 400, 1)),
                                          type = "response")) %>%
  mutate(co2 = X1,inoc = X1, co2.inoc = X1) %>%
  full_join(vcmax.beta.full) %>%
  dplyr::select(beta, co2, inoc, co2.inoc, everything(), -X1) %>%
  mutate(linetype = ifelse(co2.inoc == "amb_no.inoc", 
                           "solid", "dashed"))

##########################################################################
## Ncost:Vcmax25 plot
##########################################################################
vcmax.ncost.plot <- ggplot(data = df, 
                           aes(x = beta, y = vcmax25)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75) +
  geom_smooth(data = subset(vcmax.beta.regline, co2.inoc == "overall"),
              aes(y = emmean), size = 1.5, se = FALSE, color = "black") +
  geom_ribbon(data = subset(vcmax.beta.regline, co2.inoc == "overall"),
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
## Leaf N:total leaf N plot
##########################################################################
nmass.totalleafN <- lmer()

nmass.per.totalleafN <- emmeans()


ggplot(data = df, 
       aes(x = pnue, y = wpn)) +
  geom_jitter(aes(fill = co2.inoc),
              size = 3, shape = 21, alpha = 0.75)



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



