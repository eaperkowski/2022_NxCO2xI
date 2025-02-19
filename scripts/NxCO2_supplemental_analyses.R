## Load libraries
library(ggplot2)
library(car)
library(emmeans)
library(ggpubr)
library(lme4)
library(tidyverse)

# Load compiled datasheet
df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv", 
               na.strings = "NA") %>%
  mutate(n.trt = as.numeric(n.trt),
         rd25.vcmax25 = rd25 / vcmax25,
         inoc = factor(inoc, levels = c("no.inoc", "inoc")),
         co2 = factor(co2, levels = c("amb", "elv")),
         nod.root.ratio = nodule.biomass / root.biomass) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05)) %>%
  tidyr::unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE) 

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
## Are week 6 Vcmax25 and Jmax values different from week 7?
######################################################################
## Week 6 Vcmax25
vcmax.week <- lmer(vcmax25_wk6 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(vcmax.week))
summary(vcmax.week)
Anova(vcmax.week)
emmeans(vcmax.week, pairwise~co2)

emmeans(vcmax.week, ~1, "n.trt", at = list(n.trt = c(0, 630)))
emmeans(vcmax.week, pairwise~inoc)


## Week 6 Jmax25
jmax.week <- lmer(jmax25_wk6 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(jmax.week))
summary(jmax.week)
Anova(jmax.week)

emmeans(jmax.week, pairwise~co2)
emmeans(jmax.week, ~1, "n.trt", at = list(n.trt = c(0, 630)))
emmeans(jmax.week, pairwise~inoc)

######################################################################
## BVR analyses
######################################################################
bvr <- lmer(bvr ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(bvr, ~co2, "n.trt"))
test(emtrends(bvr, ~inoc, "n.trt"))

bvr.coefs <- data.frame(summary(bvr)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef) %>%
  dplyr::mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc", "co2:n.trt",
                              "inoc:n.trt", "co2:inoc:n.trt")) %>%
  mutate(coef = ifelse(coef <0.001 & coef >= 0,
                             "<0.001", coef)) %>%
  print(., row.names = FALSE)

bvr.table <- data.frame(Anova(bvr)) %>%
  mutate(treatment = row.names(.),
         Chisq = round(Chisq, 3),
         P_value = ifelse(Pr..Chisq. < 0.001, "<0.001",
                                round(Pr..Chisq., 3))) %>%
  full_join(bvr.coefs) %>%
  mutate(treatment = factor(
    treatment, levels = c("(Intercept)", "co2", "inoc", "n.trt", 
                          "co2:inoc", "co2:n.trt",
                          "inoc:n.trt", "co2:inoc:n.trt"))) %>%
  dplyr::select(treatment, df = Df, coef, Chisq, P_value) %>%
  arrange(treatment) %>%
  replace(is.na(.), "-")

write.csv(bvr.table, "../working_drafts/tables/NxCO2xI_tableS3_bvr.csv", row.names = FALSE)

######################################################################
## BVR plot
######################################################################
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

png("../working_drafts/figs/NxCO2xI_FigS2_bvr.png",
    width = 8, height = 4.5, units = 'in', res = 600)
bvr.plot
dev.off()

######################################################################
## Evidence against rolling rest of chlorophyll leaves
######################################################################
narea.comp <- ggplot(data = df, aes(x = narea, y = narea.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc), alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste0("atop(", ..eq.label.., ",", ..rr.label.., ")")),
               parse = TRUE, coef.digits = 4) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("inoculated",
                                "uninoculated")) +
  scale_fill_manual(values = c("#004488", "#BB5566"),
                    labels = c("ambient", "elevated")) +
  scale_x_continuous(limits = c(0,3), breaks = seq(0, 3, 1)) +
  scale_y_continuous(limits = c(0,3), breaks = seq(0, 3, 1)) +
  labs(x = expression("N"["area"]*" (gN m"^"-2"*")"),
       y = expression("N"["area_chl"]*" (gN m"^"-2"*")"),
       fill = expression("CO"[2]),
       shape = "Inoculation") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw(base_size = 20)

marea.comp <- ggplot(data = subset(df, marea.chl < 100), 
       aes(x = marea, y = marea.chl)) +
  geom_point(size = 4, aes(fill = co2, shape =  inoc), alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste0("atop(", ..eq.label.., ",", ..rr.label.., ")")),
               parse = TRUE, coef.digits = 4) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("inoculated",
                                "unnoculated")) +
  scale_fill_manual(values = c("#004488", "#BB5566"),
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
  geom_point(size = 4, aes(fill = co2, shape =  inoc), alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, size = 1) +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste0("atop(", ..eq.label.., ",", ..rr.label.., ")")),
               parse = TRUE, coef.digits = 3) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("inoculated",
                                "uninoculated")) +
  scale_fill_manual(values = c("#004488", "#BB5566"),
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
    width = 12, height = 8, units = 'in', res = 600)
ggarrange(narea.comp, nmass.comp, marea.comp,
          common.legend = TRUE, legend = "right",
          align = "hv", labels = c("(a)", "(b)", "(c)"))
dev.off()

