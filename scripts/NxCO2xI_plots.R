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
  mutate(n.trt = as.numeric(n.trt),
         rd25.vcmax25 = rd25 / vcmax25,
         inoc = factor(inoc, levels = c("no.inoc", "inoc")),
         co2 = factor(co2, levels = c("amb", "elv")),
         nod.root.ratio = nodule.biomass / root.biomass) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05)) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE) 

## Add colorblind friendly palette
nfix.cols <- c("#DDAA33", "#555555")
co2.cols <- c("#004488", "#BB5566")

full.cols <- c("#DDAA33", "#004488", "#BB5566", "#555555")

## Create blank plot as spacer plot
blank.plot <- ggplot() + 
  theme_bw() +
  theme(panel.background = element_rect(color = "white",
                                        fill = "white"),
        panel.border = element_rect(color = "white"))


##########################################################################
## Narea regression line prep
##########################################################################
narea <- lmer(narea ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(narea, ~inoc*co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
narea.co2.fert <- data.frame(emmeans(narea, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

narea.inoc.fert <- data.frame(emmeans(narea, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response"))

##########################################################################
## Narea plot
##########################################################################
narea.co2.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = narea,    
                             fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = narea.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = narea.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 3.24), breaks = seq(0, 3.2, 0.8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
narea.co2.plot


narea.fert.inoc.plot <- ggplot(data = df,
                                aes(x = n.trt,
                                    y = narea,
                                    fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = narea.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = narea.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_y_continuous(limits = c(0, 3.24), breaks = seq(0, 3.2, 0.8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
narea.fert.inoc.plot


##########################################################################
## Nmass regression line prep
##########################################################################
nmass <- lmer(nmass.focal ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nmass, ~inoc, "n.trt"))
test(emtrends(nmass, ~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nmass.co2.fert <- data.frame(emmeans(nmass, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

nmass.inoc.fert <- data.frame(emmeans(nmass, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  arrange(factor(inoc, levels = c("inoc", "no.inoc")))

##########################################################################
## Nmass plot
##########################################################################
nmass.co2.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = nmass.focal,    
                             fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = nmass.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = nmass.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["mass"]*" (gN g"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
nmass.co2.plot

nmass.fert.inoc.plot <- ggplot(data = df,
                               aes(x = n.trt,
                                   y = nmass.focal,
                                   fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = nmass.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = nmass.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["mass"]*" (gN g"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
nmass.fert.inoc.plot

##########################################################################
## Marea regression line prep
##########################################################################
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(marea, ~inoc, "n.trt"))
test(emtrends(marea, ~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
marea.co2.fert <- data.frame(emmeans(marea, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(linetype = ifelse(co2 == "elv", "solid", "dashed"))

marea.inoc.fert <- data.frame(emmeans(marea, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) 


##########################################################################
## Marea plot
##########################################################################
marea.co2.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = marea,    
                             fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = marea.co2.fert,
              aes(color = co2, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = marea.co2.fert,
              aes(fill = co2, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 15)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
marea.co2.plot


marea.fert.inoc.plot <- ggplot(data = df,
                               aes(x = n.trt,
                                   y = marea,
                                   fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = marea.inoc.fert,
              aes(color = inoc, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = marea.inoc.fert,
              aes(fill = inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 15)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
marea.fert.inoc.plot

##########################################################################
## Chlarea regression line prep
##########################################################################
chlarea <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(chlarea, ~inoc, "n.trt"))
test(emtrends(chlarea, ~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
chlarea.co2.fert <- data.frame(emmeans(chlarea, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

chlarea.inoc.fert <- data.frame(emmeans(chlarea, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response"))

##########################################################################
## Chl area plot
##########################################################################
chlarea.co2.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = chl.mmolm2,    
                             fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = chlarea.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = chlarea.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("Chl")["area"]*" (mmol m"^"-2"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
chlarea.co2.plot


chlarea.fert.inoc.plot <- ggplot(data = df,
                               aes(x = n.trt,
                                   y = chl.mmolm2,
                                   fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = chlarea.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = chlarea.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("Chl")["area"]*" (mmol m"^"-2"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25))
chlarea.fert.inoc.plot

##########################################################################
## Vcmax regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
vcmax25 <- lmer(vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax25, ~co2, "n.trt"))
test(emtrends(vcmax25, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax25.co2.fert <- data.frame(emmeans(vcmax25, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

vcmax25.inoc.fert <- data.frame(emmeans(vcmax25, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(linetype = ifelse(inoc == "inoc", "dashed", "solid"))

##########################################################################
## Vcmax plot
##########################################################################
vcmax25.co2.plot <- ggplot(data = df, 
                           aes(x = n.trt, 
                               y = vcmax25,    
                               fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = vcmax25.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = vcmax25.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
vcmax25.co2.plot


vcmax25.fert.inoc.plot <- ggplot(data = df,
                                 aes(x = n.trt,
                                     y = vcmax25,
                                     fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = vcmax25.inoc.fert,
              aes(linetype = linetype,
                  color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = vcmax25.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
vcmax25.fert.inoc.plot

##########################################################################
## Jmax regression line prep
##########################################################################
jmax25 <- lmer(jmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

test(emtrends(vcmax25, ~co2, "n.trt"))
test(emtrends(vcmax25, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jmax25.co2.fert <- data.frame(emmeans(jmax25, ~co2, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response"))

jmax25.inoc.fert <- data.frame(emmeans(jmax25, ~inoc, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response")) %>%
  mutate(linetype = ifelse(inoc == "inoc", "dashed", "solid"))


##########################################################################
## Jmax plot
##########################################################################
jmax25.co2.plot <- ggplot(data = df, 
                           aes(x = n.trt, 
                               y = jmax25,    
                               fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = jmax25.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = jmax25.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
jmax25.co2.plot


jmax25.fert.inoc.plot <- ggplot(data = df,
                                 aes(x = n.trt,
                                     y = jmax25,
                                     fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = jmax25.inoc.fert,
              aes(linetype = linetype,
                  color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = jmax25.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
jmax25.fert.inoc.plot

##########################################################################
## Jmax25:Vcmax25 regression line prep
##########################################################################
df$jmax25.vcmax25[100] <- NA
jvmax25 <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(jvmax25, ~co2, "n.trt"))
test(emtrends(jvmax25, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jvmax25.co2.fert <- data.frame(emmeans(jvmax25, ~co2, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response"))

jvmax25.inoc.fert <- data.frame(emmeans(jvmax25, ~inoc, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(linetype = ifelse(inoc == "inoc", "dashed", "solid"))


##########################################################################
## Jmax plot
##########################################################################
jvmax25.co2.plot <- ggplot(data = df, 
                          aes(x = n.trt, 
                              y = jmax25.vcmax25,    
                              fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = jvmax25.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = jvmax25.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(1.4, 2.2), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"])),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
jvmax25.co2.plot


jvmax25.fert.inoc.plot <- ggplot(data = df,
                                aes(x = n.trt,
                                    y = jmax25.vcmax25,
                                    fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = jvmax25.inoc.fert,
              aes(linetype = linetype,
                  color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = jvmax25.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(1.4, 2.2), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"])),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
jvmax25.fert.inoc.plot

##########################################################################
## Rd25 regression line prep
##########################################################################
df$rd25[df$rd25 < 0] <- NA
df$rd25[c(29, 34, 56)] <- NA

rd25 <- lmer(rd25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(rd25, ~co2, "n.trt"))
test(emtrends(rd25, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
rd25.co2.fert <- data.frame(emmeans(rd25, ~co2, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(linetype = ifelse(co2 == "elv", "dashed", "solid"))

rd25.inoc.fert <- data.frame(emmeans(rd25, ~inoc, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response")) %>%
  mutate(linetype = ifelse(inoc == "inoc", "dashed", "solid"))

##########################################################################
## Rd25 plot
##########################################################################
rd25.co2.plot <- ggplot(data = df, 
                          aes(x = n.trt, 
                              y = rd25,    
                              fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = rd25.co2.fert,
              aes(linetype = linetype, color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = rd25.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("R")["d25"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
rd25.co2.plot


rd25.fert.inoc.plot <- ggplot(data = df,
                                aes(x = n.trt,
                                    y = rd25,
                                    fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = rd25.inoc.fert,
              aes(linetype = linetype,
                  color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = rd25.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1.5)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("R")["d25"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
rd25.fert.inoc.plot

##########################################################################
## Prop leaf N in photosynthesis  regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$p.photo[df$p.photo > 1] <- NA
p.photo <- lmer(p.photo ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(p.photo, ~inoc, "n.trt"))
test(emtrends(p.photo, pairwise~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.photo.co2.fert <- data.frame(emmeans(p.photo, ~co2, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(linetype = ifelse(co2 == "amb", "dashed", "solid"))

p.photo.inoc.fert <- data.frame(emmeans(p.photo, ~inoc, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response")) %>%
  mutate(linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))

##########################################################################
## Prop leaf N in photosynthesis plot
##########################################################################
p.photo.co2.plot <- ggplot(data = df, 
                        aes(x = n.trt, 
                            y = p.photo,    
                            fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = p.photo.co2.fert,
              aes(linetype = linetype, color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = p.photo.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["photo"]*" (g g"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
p.photo.co2.plot


p.photo.fert.inoc.plot <- ggplot(data = df,
                              aes(x = n.trt,
                                  y = p.photo,
                                  fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = p.photo.inoc.fert,
              aes(linetype = linetype,
                  color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = p.photo.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["photo"]*" (g g"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
p.photo.fert.inoc.plot

##########################################################################
## Prop leaf N in structure regression line prep
##########################################################################
p.str <- lmer(log(p.structure) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(p.str, ~co2, "n.trt"))
test(emtrends(p.str, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.str.co2.fert <- data.frame(emmeans(p.str, ~co2, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response"))

p.str.inoc.fert <- data.frame(emmeans(p.str, ~inoc, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response"))

##########################################################################
## Prop leaf N in structure plot
##########################################################################
p.str.co2.plot <- ggplot(data = df, 
                           aes(x = n.trt, 
                               y = p.structure,    
                               fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = p.str.co2.fert,
              aes(color = co2, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = p.str.co2.fert,
              aes(fill = co2, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["structure"]*" (g g"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
p.str.co2.plot


p.str.fert.inoc.plot <- ggplot(data = df,
                                 aes(x = n.trt,
                                     y = p.structure,
                                     fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = p.str.inoc.fert,
              aes(color = inoc, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = p.str.inoc.fert,
              aes(fill = inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["structure"]*" (g g"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
p.str.fert.inoc.plot

##########################################################################
## PNUE regression line prep
##########################################################################
df$pnue[41] <- NA
pnue <- lmer(pnue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(pnue, ~co2, "n.trt"))
test(emtrends(pnue, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
pnue.co2.fert <- data.frame(emmeans(pnue, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

pnue.inoc.fert <- data.frame(emmeans(pnue, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))


##########################################################################
## PNUE plot
##########################################################################
pnue.co2.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = pnue,    
                             fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = pnue.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = pnue.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("PNUE ("*mu*"mol CO"["2"]*" mol N"^"-1"*" s"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
pnue.co2.plot


pnue.fert.inoc.plot <- ggplot(data = df,
                               aes(x = n.trt,
                                   y = pnue,
                                   fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = pnue.inoc.fert,
              aes(linetype = linetype, color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = pnue.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("PNUE ("*mu*"mol CO"["2"]*" mol N"^"-1"*" s"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
pnue.fert.inoc.plot

##########################################################################
## iWUE regression line prep - to be replaced with chi?
##########################################################################
iwue <- lmer(iwue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(iwue, ~co2, "n.trt"))
test(emtrends(iwue, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
iwue.co2.fert <- data.frame(emmeans(iwue, ~co2, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response"))

iwue.inoc.fert <- data.frame(emmeans(iwue, ~inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

##########################################################################
## iWUE plot
##########################################################################
iwue.co2.plot <- ggplot(data = df, 
                        aes(x = n.trt, 
                            y = iwue,    
                            fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = iwue.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = iwue.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(25, 125), breaks = seq(25, 125, 25)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("iWUE ("*mu*"mol CO"["2"]*" mol"^"-1"*" H"["2"]*"O)")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
iwue.co2.plot


iwue.fert.inoc.plot <- ggplot(data = df,
                              aes(x = n.trt,
                                  y = iwue,
                                  fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = iwue.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = iwue.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(25, 125), breaks = seq(25, 125, 25)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("iWUE ("*mu*"mol CO"["2"]*" mol"^"-1"*" H"["2"]*"O)")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
iwue.fert.inoc.plot

##########################################################################
## Narea:gs regression line prep
##########################################################################
narea.gs <- lmer(log(narea.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(narea.gs, ~co2, "n.trt"))
test(emtrends(narea.gs, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
narea.gs.co2.fert <- data.frame(emmeans(narea.gs, ~co2, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response"))

narea.gs.inoc.fert <- data.frame(emmeans(narea.gs, ~inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))

##########################################################################
## Narea:gs plot
##########################################################################
narea.gs.co2.plot <- ggplot(data = df, 
                            aes(x = n.trt, 
                                y = narea.gs,    
                                fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = narea.gs.co2.fert,
              aes(color = co2, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = narea.gs.co2.fert,
              aes(fill = co2, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 32, 8)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*": "*italic("g")["s"]*" (gN s mol"^"-1"*" H"["2"]*"O)")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
narea.gs.co2.plot


narea.gs.fert.inoc.plot <- ggplot(data = df,
                              aes(x = n.trt,
                                  y = narea.gs,
                                  fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = narea.gs.inoc.fert,
              aes(color = inoc, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = narea.gs.inoc.fert,
              aes(fill = inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 32, 8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*": "*italic("g")["s"]*" (gN s mol"^"-1"*" H"["2"]*"O)")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
narea.gs.fert.inoc.plot

##########################################################################
## Vcmax:gs regression line prep
##########################################################################
vcmax.gs <- lmer(log(vcmax.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax.gs, ~co2, "n.trt"))
test(emtrends(vcmax.gs, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax.gs.co2.fert <- data.frame(emmeans(vcmax.gs, ~co2, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response"))

vcmax.gs.inoc.fert <- data.frame(emmeans(vcmax.gs, ~inoc, "n.trt",
                                         at = list(n.trt = seq(0, 630, 5)),
                                         type = "response")) %>%
  mutate(linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))

##########################################################################
## Vcmax:gs plot
##########################################################################
vcmax.gs.co2.plot <- ggplot(data = df, 
                            aes(x = n.trt, 
                                y = vcmax.gs,    
                                fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = vcmax.gs.co2.fert,
              aes(color = co2, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = vcmax.gs.co2.fert,
              aes(fill = co2, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 200)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*": "*italic("g")["s"]*" ("*mu*"mol CO"["2"]*" mol"^"-1"*" H"["2"]*"O)")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
vcmax.gs.co2.plot


vcmax.gs.fert.inoc.plot <- ggplot(data = df,
                                  aes(x = n.trt,
                                      y = vcmax.gs,
                                      fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = vcmax.gs.inoc.fert,
              aes(color = inoc, y = response), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = vcmax.gs.inoc.fert,
              aes(fill = inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 200)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*": "*italic("g")["s"]*" ("*mu*"mol CO"["2"]*" mol"^"-1"*" H"["2"]*"O)")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        axis.title.y = element_text(size = 14)) +
  guides(linetype = "none")
vcmax.gs.fert.inoc.plot


##########################################################################
## Total leaf area regression line prep
##########################################################################
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(tla, ~co2, "n.trt"))
test(emtrends(tla, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
tla.co2.fert <- data.frame(emmeans(tla, ~co2, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response"))

tla.inoc.fert <- data.frame(emmeans(tla, ~inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response"))


##########################################################################
## Total leaf area plot
##########################################################################
tla.co2.plot <- ggplot(data = df, 
                            aes(x = n.trt, 
                                y = tla,    
                                fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = tla.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = tla.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 300)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
tla.co2.plot


tla.fert.inoc.plot <- ggplot(data = df,
                                  aes(x = n.trt,
                                      y = tla,
                                      fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = tla.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = tla.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 300)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
tla.fert.inoc.plot

##########################################################################
## Total biomass regression line prep
##########################################################################
tbio <- lmer(total.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(tbio, ~co2, "n.trt"))
test(emtrends(tbio, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
tbio.co2.fert <- data.frame(emmeans(tbio, ~co2, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response"))

tbio.inoc.fert <- data.frame(emmeans(tbio, ~inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response"))

##########################################################################
## Total biomass plot
##########################################################################
tbio.co2.plot <- ggplot(data = df, 
                       aes(x = n.trt, 
                           y = total.biomass,    
                           fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = tbio.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = tbio.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Total biomass (g)",
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
tbio.co2.plot


tbio.fert.inoc.plot <- ggplot(data = df,
                             aes(x = n.trt,
                                 y = total.biomass,
                                 fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = tbio.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = tbio.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Total biomass (g)",
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
tbio.fert.inoc.plot

##########################################################################
## Ncost regression line prep
##########################################################################
df$ncost[c(100, 101)] <- NA
df$ncost[c(38, 103)] <- NA
df$ncost[32] <- NA

ncost <- lmer(ncost ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(ncost, ~co2, "n.trt"))
test(emtrends(ncost, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
ncost.co2.fert <- data.frame(emmeans(ncost, ~co2, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response"))

ncost.inoc.fert <- data.frame(emmeans(ncost, ~inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

##########################################################################
## Ncost plot
##########################################################################
ncost.co2.plot <- ggplot(data = df, 
                        aes(x = n.trt, 
                            y = ncost,    
                            fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = ncost.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = ncost.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
ncost.co2.plot


ncost.fert.inoc.plot <- ggplot(data = df,
                              aes(x = n.trt,
                                  y = ncost,
                                  fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = ncost.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = ncost.inoc.fert,
              aes(fill = inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = nfix.cols,
                    labels = c("Uninoculated",
                               "Inoculated")) +
  scale_color_manual(values = nfix.cols,
                     labels = c("Uninoculated",
                                "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
ncost.fert.inoc.plot

##########################################################################
## Root nodule biomass regression line prep
##########################################################################
df$nodule.biomass[df$nod.root.ratio > 0.05 & df$inoc == "no.inoc"] <- NA
df$nodule.biomass[80] <- NA

nod <- lmer(sqrt(nodule.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nod, ~co2, "n.trt"))
test(emtrends(nod, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nod.regline <- data.frame(emmeans(nod, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))

##########################################################################
## Root nodule biomass regression line prep
##########################################################################
nod.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = nodule.biomass,    
                             fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = nod.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = nod.regline,
              aes(fill = co2.inoc, y = response,
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = full.cols,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = full.cols,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(-0.01, 0.6), breaks = seq(0, 0.6, 0.15)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Nodule biomass (g)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
nod.plot

##########################################################################
## Root nodule:root biomass regression line prep
##########################################################################
df$nod.root.ratio <- df$nodule.biomass / df$root.biomass
df$nod.root.ratio[df$nod.root.ratio > 0.05 & df$inoc == "no.inoc"] <- NA

nodroot <- lmer(sqrt(nod.root.ratio) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nodroot, ~co2*inoc, "n.trt"))
test(emtrends(nodroot, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nodroot.regline <- data.frame(emmeans(nodroot, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))

##########################################################################
## Root nodule:root biomass regression line prep
##########################################################################
nodroot.plot <- ggplot(data = df, 
                       aes(x = n.trt, 
                           y = nod.root.ratio,    
                           fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = nodroot.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = nodroot.regline,
              aes(fill = co2.inoc, y = response,
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = full.cols,
                    labels = c("Ambient, inoculated",
                               "Ambient, uninoculated",
                               "Elevated, inoculated",
                               "Elevated, uninoculated")) +
  scale_color_manual(values = full.cols,
                     labels = c("Ambient, inoculated",
                                "Ambient, uninoculated",
                                "Elevated, inoculated",
                                "Elevated, uninoculated")) +
  scale_y_continuous(limits = c(-0.0001, 0.4), breaks = seq(0, 0.4, 0.1)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Nodule: root biomass",
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
nodroot.plot

##########################################################################
## %Ndfa
##########################################################################
fake <- data.frame(x = seq(0,630,1), y = seq(0,630,1))
label <- data.frame(x = 315, y = 0.5, text = "placeholder for %Ndfa figure")

ndfa.plot <- ggplot(data = fake, aes(x=x,y=y)) +
  geom_text(data = label, aes(label=text), fontface = "bold", size = 5) +
  scale_x_continuous(limits = c(0, 650), breaks = seq(0,600,200)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
  labs(x = expression(bold("Soil N fertilization (ppm)")),
       y = expression(bold("%N"["dfa"]))) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1.25))

##########################################################################
## Figure 1: leaf N plots
##########################################################################
legend_co2_plots <- get_legend(narea.co2.plot)
legend_fert_plots <- get_legend(narea.fert.inoc.plot)
legend_blank_plot <- get_legend(blank.plot)

png("../working_drafts/figs/NxCO2xI_fig1_leafN.png",
    height = 16, width = 12, units = "in", res = 600)
ggarrange(ggarrange(narea.co2.plot, nmass.co2.plot, marea.co2.plot, 
                    chlarea.co2.plot, ncol = 1, nrow = 4, 
                    align = "hv", legend = "none", 
                    labels = c("A", "C", "E", "G"), font.label = list(size = 18)),
          ggarrange(narea.fert.inoc.plot, nmass.fert.inoc.plot, 
                    marea.fert.inoc.plot, chlarea.fert.inoc.plot, ncol = 1, 
                    nrow = 4, align = "hv", legend = "none",
                    labels = c("B", "D", "F", "H"), font.label = list(size = 18)), 
          ggarrange(legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot, legend_blank_plot,
                    legend_co2_plots, legend_fert_plots, 
                    legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot, legend_blank_plot,
                    nrow = 12, align = "hv"),
          ncol = 3, align = "hv", legend = "none", widths = c(0.4, 0.4, 0.2))
dev.off()

##########################################################################
## Figure 2: leaf physiology plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig2_photo.png",
    height = 12, width = 12, units = "in", res = 600)
ggarrange(ggarrange(vcmax25.co2.plot, jmax25.co2.plot,
                    jvmax25.co2.plot, ncol = 1, nrow = 3, 
                    align = "hv", legend = "none", 
                    labels = c("A", "C", "E", "G"), font.label = list(size = 18)),
          ggarrange(vcmax25.fert.inoc.plot, jmax25.fert.inoc.plot, 
                    jvmax25.fert.inoc.plot, ncol = 1, 
                    nrow = 3, align = "hv", legend = "none",
                    labels = c("B", "D", "F", "H"), font.label = list(size = 18)), 
          ggarrange(legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot, 
                    legend_co2_plots, legend_fert_plots, 
                    legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot,
                    nrow = 10, align = "hv"),
          ncol = 3, align = "hv", legend = "none", widths = c(0.4, 0.4, 0.2))

dev.off()

##########################################################################
## Figure 3: propN photosynthesis/structure
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig3_propN.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(ggarrange(p.photo.co2.plot, p.str.co2.plot, ncol = 1, nrow = 2, 
                    align = "hv", legend = "none", 
                    labels = c("A", "C"), font.label = list(size = 18)),
          ggarrange(p.photo.fert.inoc.plot, p.str.fert.inoc.plot, ncol = 1, 
                    nrow = 2, align = "hv", legend = "none",
                    labels = c("B", "D"), font.label = list(size = 18)), 
          ggarrange(legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_co2_plots, legend_fert_plots, 
                    legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    nrow = 8, align = "hv"),
          ncol = 3, align = "hv", legend = "none", widths = c(0.4, 0.4, 0.2))
dev.off()

##########################################################################
## Figure 4: PNUE/iWUE tradeoffs
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig4_PNUE_iWUE.png",
    height = 16, width = 12, units = "in", res = 600)
ggarrange(ggarrange(pnue.co2.plot, iwue.co2.plot,
                    narea.gs.co2.plot, vcmax.gs.co2.plot, ncol = 1, nrow = 4, 
                    align = "hv", legend = "none", 
                    labels = c("A", "C", "E", "G"), font.label = list(size = 18)),
          ggarrange(pnue.fert.inoc.plot, iwue.fert.inoc.plot, 
                    narea.gs.fert.inoc.plot, vcmax.gs.fert.inoc.plot,
                    ncol = 1, nrow = 4, align = "hv", legend = "none",
                    labels = c("B", "D", "F", "H"), font.label = list(size = 18)), 
          ggarrange(legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot, legend_blank_plot,
                    legend_co2_plots, legend_fert_plots, 
                    legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot, legend_blank_plot,
                    nrow = 12, align = "hv"),
          ncol = 3, align = "hv", legend = "none", widths = c(0.4, 0.4, 0.2))
dev.off()

##########################################################################
## Figure 5: whole plant plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig5_wholePlant.png",
    height = 12, width = 12, units = "in", res = 600)
ggarrange(ggarrange(tla.co2.plot, tbio.co2.plot,
                    ncost.co2.plot, ncol = 1, nrow = 3, 
                    align = "hv", legend = "none", 
                    labels = c("A", "C", "E"), font.label = list(size = 18)),
          ggarrange(tla.fert.inoc.plot, tbio.fert.inoc.plot, 
                    ncost.fert.inoc.plot, ncol = 1, 
                    nrow = 3, align = "hv", legend = "none",
                    labels = c("B", "D", "F"), font.label = list(size = 18)), 
          ggarrange(legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot, 
                    legend_co2_plots, legend_fert_plots, 
                    legend_blank_plot, legend_blank_plot, legend_blank_plot,
                    legend_blank_plot,
                    nrow = 10, align = "hv"),
          ncol = 3, align = "hv", legend = "none", widths = c(0.4, 0.4, 0.2))
dev.off()

##########################################################################
## Figure 6: nitrogen fixation plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig6_nFix.png",
    height = 12, width = 7, units = "in", res = 600)
ggarrange(nod.plot, nodroot.plot, ndfa.plot, 
          ncol = 1, nrow = 3, align = "hv", common.legend = TRUE,
          legend = "right", labels = "AUTO", font.label = list(size = 18))
dev.off()





