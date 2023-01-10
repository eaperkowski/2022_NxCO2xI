## Load libraries
library(ggplot2)
library(car)
library(emmeans)
library(ggpubr)

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


######################################################################
## Belowground carbon biomass
######################################################################
cbg <- lmer(cbg ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(cbg, ~co2, "n.trt"))
test(emtrends(cbg, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
cbg.co2.fert <- data.frame(emmeans(cbg, ~co2, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response"))

cbg.inoc.fert <- data.frame(emmeans(cbg, ~inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response"))

##########################################################################
## cbg plot
##########################################################################
cbg.co2.plot <- ggplot(data = df, 
                         aes(x = n.trt, 
                             y = cbg,    
                             fill = co2)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = cbg.co2.fert,
              aes(color = co2, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = cbg.co2.fert,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              size = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient",
                               "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient",
                                "Elevated")) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("C")["bg"]*" (gC)")),
       fill = expression(bold("CO"["2"])), color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
cbg.co2.plot


cbg.fert.inoc.plot <- ggplot(data = df,
                               aes(x = n.trt,
                                   y = cbg,
                                   fill = inoc)) +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data = cbg.inoc.fert,
              aes(color = inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = cbg.inoc.fert,
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
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("C")["bg"]*" (gC)")),
       fill = "Inoculation", color = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
cbg.fert.inoc.plot






######################################################################
## BVR
######################################################################
bvr <- lmer(bvr ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nodroot, ~co2, "n.trt"))
test(emtrends(nodroot, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
bvr.regline <- data.frame(emmeans(bvr, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc))



bvr.plot <- ggplot(data = df, 
                   aes(x = n.trt, 
                       y = bvr,    
                       fill = co2.inoc)) +
  geom_hline(yintercept = 1, col = "black", linetype = "dotted") +
  geom_hline(yintercept = 2, col = "black", linetype = "dashed") +
  geom_jitter(size = 3, alpha = 0.75, shape = 21) +
  geom_smooth(data =bvr.regline,
              aes(color = co2.inoc, y = emmean), 
              size = 1.5, se = FALSE) +
  geom_ribbon(data = bvr.regline,
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
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("BVR (g L"^"-1"*")")),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25)) +
  guides(linetype = "none")
bvr.plot

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

