## Load libraries
library(ggplot2)
library(car)
library(emmeans)
library(ggpubr)

## Load file
df <- read.csv("../data_sheets/NxCO2_datasheet.csv")

## Preliminary models for LEMONTREE meeting
df$vcmax25[91] <- NA

vcmax25 <- lm(vcmax25 ~ as.factor(co2) * inoc * as.numeric(n.trt),
              data = df)
summary(vcmax25)
Anova(vcmax25)

shapiro.test(residuals(vcmax25))
outlierTest(vcmax25)

## Preliminary plots for LEMONTREE meeting
inoc.labs <- c("Not inoculated", "Inoculated")
names(inoc.labs) <- c("n", "y")


vcmax25 <- ggplot(data = subset(df, !is.na(inoc)), 
       aes(x = as.numeric(n.trt), y = vcmax25)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 30)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  labs(x = NULL,
       y = expression(italic("V")["cmax25"]*" (μmol m"^"-2"*"s"^"-1"*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 22)

jmax25 <- ggplot(data = subset(df, !is.na(inoc)), 
       aes(x = as.numeric(n.trt), y = jmax25)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("J")["max25"]*" (μmol m"^"-2"*"s"^"-1"*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  theme_bw(base_size = 22)


jv.rat <- ggplot(data = subset(df, !is.na(inoc)), 
       aes(x = as.numeric(n.trt), y = jmax25.vcmax25)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  scale_y_continuous(limits = c(1.4,2), breaks = seq(1.4,2,0.2)) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression(italic("J")["max25"]*":"*italic("V")["cmax25"]),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 22)

tla <- ggplot(data = subset(df, !is.na(inoc)), 
              aes(x = as.numeric(n.trt), y = tla)) +
  geom_point(shape = 21, size = 5, alpha = 0.9, aes(fill = co2)) +
  geom_smooth(method = "lm", aes(color = co2)) +
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 300)) +
  labs(x = "Soil nitrogen fertilization (ppm N)",
       y = expression("Total leaf area (cm"^2*")"),
       fill = expression("CO"["2"]*" treatment")) +
  facet_grid(~inoc, labeller = labeller(inoc = inoc.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 22)

png("../working_drafts/NxCO2_LEMONTREE_fig_vcmax_jvrat.png",
    width = 12, height = 12, units = 'in', res = 600)
ggarrange(vcmax25, jv.rat, nrow = 2, align = "hv", common.legend = TRUE,
          legend = "right")
dev.off()

png("../working_drafts/NxCO2_LEMONTREE_fig_tla.png",
    width = 12, height = 6, units = 'in', res = 600)
tla
dev.off()

