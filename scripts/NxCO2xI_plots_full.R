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
narea.regline <- data.frame(emmeans(narea, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc))

##########################################################################
## Narea plot
##########################################################################
narea.plot <- ggplot(data = df, 
                     aes(x = n.trt, 
                         y = narea,    
                         fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = narea.regline,
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = narea.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 3.24), breaks = seq(0, 3.2, 0.8)) +
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
nmass.regline <- data.frame(emmeans(nmass, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"))

##########################################################################
## Nmass plot
##########################################################################
nmass.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = nmass.focal,    
                             fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = nmass.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = nmass.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["mass"]*" (gN g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
nmass.plot


##########################################################################
## Marea regression line prep
##########################################################################
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(marea, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
marea.regline <- data.frame(emmeans(marea, ~co2*inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"))

##########################################################################
## Marea plot
##########################################################################
marea.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = marea,    
                             fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = marea.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = marea.regline,
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
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 15)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
marea.plot

##########################################################################
## Chlarea regression line prep
##########################################################################
chlarea <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(chlarea, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
chlarea.regline <- data.frame(emmeans(chlarea, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc))

##########################################################################
## Chl area plot
##########################################################################
chlarea.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = chl.mmolm2,    
                             fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = chlarea.regline,
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = chlarea.regline,
              aes(fill = co2.inoc, y = emmean, 
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
chlarea.plot


##########################################################################
## Vcmax regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
vcmax25 <- lmer(vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax25.regline <- data.frame(emmeans(vcmax25, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "inoc", "dashed", "solid"))

##########################################################################
## Vcmax plot
##########################################################################
vcmax25.plot <- ggplot(data = df, 
                           aes(x = n.trt, 
                               y = vcmax25,    
                               fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = vcmax25.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = vcmax25.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
vcmax25.plot

##########################################################################
## Jmax regression line prep
##########################################################################
jmax25 <- lmer(jmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jmax25.regline <- data.frame(emmeans(jmax25, ~co2*inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "inoc", "dashed", "solid"))


##########################################################################
## Jmax plot
##########################################################################
jmax25.plot <- ggplot(data = df, 
                           aes(x = n.trt, 
                               y = jmax25,    
                               fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = jmax25.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = jmax25.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
jmax25.plot



##########################################################################
## Jmax25:Vcmax25 regression line prep
##########################################################################
df$jmax25.vcmax25[100] <- NA
jvmax25 <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(jvmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jvmax25.regline <- data.frame(emmeans(jvmax25, ~co2*inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "inoc", "dashed", "solid"))


##########################################################################
## Jmax plot
##########################################################################
jvmax25.plot <- ggplot(data = df, 
                          aes(x = n.trt, 
                              y = jmax25.vcmax25,    
                              fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = jvmax25.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = jvmax25.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(1.4, 2.2), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"])),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
jvmax25.plot

##########################################################################
## Rd25 regression line prep
##########################################################################
df$rd25[df$rd25 < 0] <- NA
df$rd25[c(29, 34, 56)] <- NA

rd25 <- lmer(rd25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(rd25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
rd25.regline <- data.frame(emmeans(rd25, ~co2*inoc, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_no.inoc", "solid", "dashed"))

##########################################################################
## Rd25 plot
##########################################################################
rd25.plot <- ggplot(data = df, 
                          aes(x = n.trt, 
                              y = rd25,    
                              fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = rd25.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = rd25.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("R")["d25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
rd25.plot

##########################################################################
## gsw regression line prep
##########################################################################
df$gsw[69] <- NA

gsw <- lmer(gsw ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(gsw, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
gsw.regline <- data.frame(emmeans(gsw, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_no.inoc" |
                             co2.inoc == "amb_inoc", "dashed", "solid"))

##########################################################################
## gsw plot
##########################################################################
gsw.plot <- ggplot(data = df, 
                    aes(x = n.trt, 
                        y = gsw,    
                        fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = gsw.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = gsw.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.15)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("g")["sw"]*" (mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
gsw.plot


##########################################################################
## gsw regression line prep
##########################################################################
l <- lmer(stomlim ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(l, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
l.regline <- data.frame(emmeans(l, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"))

##########################################################################
## gsw plot
##########################################################################
l.plot <- ggplot(data = df, 
                   aes(x = n.trt, 
                       y = stomlim,    
                       fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = l.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = l.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.15)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Stomatal limitation")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
l.plot


##########################################################################
## Prop leaf N in photosynthesis  regression line prep
##########################################################################
## Copy removed outliers and lmer fxn
df$p.photo[df$p.photo > 1] <- NA
p.photo <- lmer(p.photo ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(p.photo, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.photo.regline <- data.frame(emmeans(p.photo, ~co2*inoc, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_no.inoc", "dashed", "solid"))

##########################################################################
## Prop leaf N in photosynthesis plot
##########################################################################
p.photo.plot <- ggplot(data = df, 
                        aes(x = n.trt, 
                            y = p.photo,    
                            fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = p.photo.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = p.photo.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["photo"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
p.photo.plot


##########################################################################
## Prop leaf N in structure regression line prep
##########################################################################
p.str <- lmer(log(p.structure) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(p.str, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
p.str.regline <- data.frame(emmeans(p.str, ~co2*inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"))

##########################################################################
## Prop leaf N in structure plot
##########################################################################
p.str.plot <- ggplot(data = df, 
                           aes(x = n.trt, 
                               y = p.structure,    
                               fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = p.str.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = p.str.regline,
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
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic(rho)["structure"]*" (g g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
p.str.plot

##########################################################################
## Total leaf area regression line prep
##########################################################################
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(tla, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
tla.regline <- data.frame(emmeans(tla, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc))


##########################################################################
## Total leaf area plot
##########################################################################
tla.plot <- ggplot(data = df, 
                   aes(x = n.trt, 
                       y = tla,
                       fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = tla.regline,
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = tla.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 300)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
tla.plot

##########################################################################
## Total biomass regression line prep
##########################################################################
tbio <- lmer(total.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(tbio, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
tbio.regline <- data.frame(emmeans(tbio, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc))


##########################################################################
## Total biomass plot
##########################################################################
tbio.plot <- ggplot(data = df, 
                       aes(x = n.trt, 
                           y = total.biomass,    
                           fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = tbio.regline,
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = tbio.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Total biomass (g)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
tbio.plot


##########################################################################
## Ncost regression line prep
##########################################################################
df$ncost[c(100, 101)] <- NA
df$ncost[c(38, 103)] <- NA
df$ncost[32] <- NA

ncost <- lmer(ncost ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(ncost, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
ncost.regline <- data.frame(emmeans(ncost, ~co2*inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"))

##########################################################################
## Ncost plot
##########################################################################
ncost.plot <- ggplot(data = df, 
                        aes(x = n.trt, 
                            y = ncost,    
                            fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = ncost.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = ncost.regline,
              aes(fill = co2.inoc, y = emmean, 
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
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
ncost.plot

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
df$ndfa[c(38, 85, 101, 103)] <- NA
ndfa <- lmer(sqrt(ndfa) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(ndfa, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
ndfa.regline <- data.frame(emmeans(ndfa, ~co2*inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "no.inoc", "dashed", "solid"))

ndfa.plot <- ggplot(data = df, 
                    aes(x = n.trt, 
                        y = ndfa,    
                        fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = ndfa.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = ndfa.regline,
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
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("%N"["dfa"])),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
ndfa.plot

##########################################################################
## Exploratory: effects of treatment combos on beta
##########################################################################
beta <- lmer(log(beta) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(beta, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
beta.regline <- data.frame(emmeans(beta, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"))

beta.plot <- ggplot(data = df, 
                    aes(x = chi, 
                        y = narea,    
                        fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = beta.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = beta.regline,
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
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 100)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(beta)),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
beta.plot

##########################################################################
## Exploratory: effects of treatment combos on beta
##########################################################################
beta <- lmer(log(beta) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(beta, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
beta.regline <- data.frame(emmeans(beta, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"))

beta.plot <- ggplot(data = df, 
                    aes(x = n.trt, 
                        y = beta,    
                        fill = co2.inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = beta.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = beta.regline,
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
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 100)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(beta)),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
beta.plot

##########################################################################
## Figure 1: leaf N plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig1_leafN.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(narea.plot, nmass.plot, marea.plot, chlarea.plot, 
          ncol = 2, nrow = 2, align = "hv", legend = "right",
          common.legend = TRUE, labels = c("(a)", "(b)", "(c)", "(d)"),
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 2: leaf physiology plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig2_photo.png",
    height = 12, width = 12, units = "in", res = 600)
ggarrange(vcmax25.plot, jmax25.plot, jvmax25.plot, 
          rd25.plot, gsw.plot, l.plot,
          ncol = 2, nrow = 3, align = "hv", legend = "right",
          common.legend = TRUE, font.label = list(size = 18), 
          labels = c("(a)", "(b)", "(c)", "(d)","(e)", "(f)"))
dev.off()

##########################################################################
## Figure 3: propN photosynthesis/structure
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig3_propN.png",
    height = 4, width = 12, units = "in", res = 600)
ggarrange(p.photo.plot, p.str.plot, ncol = 2, nrow = 1,
          align = "hv", legend = "right", common.legend = TRUE,
          labels = c("(a)", "(b)"), font.label = list(size = 18))
dev.off()

##########################################################################
## Figure 4: whole plant plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_fig5_wholePlant.png",
    height = 8, width = 12, units = "in", res = 600)
ggarrange(tla.plot, tbio.plot, ncost.plot, ndfa.plot,
          ncol = 2, nrow = 2, align = "hv", legend = "right",
          labels = c("(a)", "(b)", "(c)", "(d)"), common.legend = TRUE,
          font.label = list(size = 18))
dev.off()

##########################################################################
## Figure SX: nitrogen fixation plots
##########################################################################
png("../working_drafts/figs/NxCO2xI_figS3_nFix.png",
    height = 4, width = 12, units = "in", res = 600)
ggarrange(nod.plot, nodroot.plot,
          align = "hv", common.legend = TRUE,
          nrow = 1, ncol = 2,
          legend = "right", labels = c("(a)", "(b)"), 
          font.label = list(size = 18))
dev.off()





