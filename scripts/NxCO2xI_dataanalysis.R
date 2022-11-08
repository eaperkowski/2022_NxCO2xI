##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Read in compiled data file
df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  mutate(n.trt = as.numeric(n.trt))

## Some light analyses while I construct actual chamber environment script
par.chamber <- data.frame(percent = seq(0, 100, 10),
                          ch1 = c(0, 107, 222, 327, 443, 552, 615, 705, 877, 1004, 1225),
                          ch2 = c(0, 106, 220, 325, 440, 547, 610, 700, 870, 996,  1215),
                          ch3 = c(0, 119, 220, 324, 439, 548, 611, 700, 871, 998,  1217),
                          ch4 = c(0, 123, 228, 337, 456, 568, 634, 727, 905, 1037, 1265),
                          ch5 = c(0, 112, 233, 344, 466, 581, 649, 744, 926, 1061, 1295),
                          ch6 = c(0, 107, 222, 328, 444, 552, 616, 706, 879, 1006, 1228))

summary(lm(ch1 ~ percent, data = par.chamber))
summary(lm(ch2 ~ percent, data = par.chamber))
summary(lm(ch3 ~ percent, data = par.chamber))
summary(lm(ch4 ~ percent, data = par.chamber))
summary(lm(ch5 ~ percent, data = par.chamber))
summary(lm(ch6 ~ percent, data = par.chamber))

# Mean and standard deviation for maximum PAR
mean(as.numeric(par.chamber[11, c(2:7)]))
sd(as.numeric(par.chamber[11, c(2:7)]))

## Daily 16:8 light availability during day
day.par <- data.frame(chamber = seq(1,6,1), 
                      day.par = c(((11.45*25 - 20.5455)*1.5 + 
                                     (11.45*50 - 20.5455)*1.5 + 
                                     (11.45*75 - 20.5455)*1.5 +
                                     (1225*11.5))/16,
                                  ((11.3682*25 -20.3182)*1.5 + 
                                     (11.3682*50 -20.3182)*1.5 + 
                                     (11.3682*75 -20.3182)*1.5 +
                                     (1215*11.5))/16,
                                  ((11.344*25 -17.454)*1.5 + 
                                     (11.344*50 -17.454)*1.5 + 
                                     (11.344*75 -17.454)*1.5 +
                                     (1217*11.5))/16,
                                  ((11.7909*25 -18.6364)*1.5 + 
                                     (11.7909*50 -18.6364)*1.5 + 
                                     (11.7909*75 -18.6364)*1.5 +
                                     (1265*11.5))/16,
                                  ((12.1209*25 -23.2273)*1.5 + 
                                     (12.1209*50 -23.2273)*1.5 + 
                                     (12.1209*75 -23.2273)*1.5 +
                                     (1295*11.5))/16,
                                  ((11.4864*25 -20.8636)*1.5 + 
                                     (11.4864*50 -20.8636)*1.5 + 
                                     (11.4864*75 -20.8636)*1.5 +
                                     (1228*11.5))/16))

mean(day.par$day.par)
sd(day.par$day.par)

##########################################################################
## Leaf nitrogen content (Narea)
##########################################################################
narea <- lmer(narea ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(narea)
qqnorm(residuals(narea))
qqline(residuals(narea))
densityPlot(residuals(narea))
shapiro.test(residuals(narea))
outlierTest(narea)

# Model results
summary(narea)
Anova(narea)
r.squaredGLMM(narea)

# Post-hoc tests
# Interaction between CO2 and n.trt
test(emtrends(narea, pairwise~co2, "n.trt")) 
# Stronger stimulation in Narea with increasing soil N under ambient CO2

# Interaction between inoc and n.trt
test(emtrends(narea, pairwise~inoc, "n.trt"))
# Stronger stimulation in Narea with increasing soil N in non-inoculated pots

# Interaction between CO2 and inoc
emmeans(narea, pairwise~inoc*co2)
# Stronger stimulation in Narea in inoculated pots grown under elevated CO2

# Individual effect of co2
emmeans(narea, pairwise~co2)
# Narea is generally greater under ambient CO2

# Individual effect of inoc
emmeans(narea, pairwise~inoc)
# Narea is generally larger in inoculated pots

# Individual effect of n.trt
emtrends(narea, ~1, "n.trt")
# Narea increases with N fertilization

##########################################################################
## Marea (LMA)
##########################################################################
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
qqnorm(residuals(marea))
qqline(residuals(marea))
densityPlot(residuals(marea))
shapiro.test(residuals(marea))
outlierTest(marea)

# Model results
summary(marea)
Anova(marea)
r.squaredGLMM(marea)

# Post-hoc tests
# Interaction between CO2 and n.trt
test(emtrends(marea, pairwise~co2, "n.trt")) 
# Stronger stimulation in Marea with increasing soil N under elevated CO2

# Interaction between inoc and n.trt
test(emtrends(marea, pairwise~inoc, "n.trt"))
# Stronger stimulation in Narea with increasing soil N in non-inoculated pots

# Individual effect of co2
emmeans(marea, pairwise~co2)
# Marea is larger in elevated CO2 treatment

# Individual effect of inoc
emmeans(marea, pairwise~inoc)
# Marea is larger in inoculated pots

# Individual effect of n.trt
emtrends(marea, ~1, "n.trt")
# Increasing N fertilization increases Marea

##########################################################################
## Nmass
##########################################################################
df$nmass.focal[c(110, 111, 114)] <- NA

nmass <- lmer(nmass.focal ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
qqnorm(residuals(nmass))
qqline(residuals(nmass))
densityPlot(residuals(nmass))
shapiro.test(residuals(nmass))
outlierTest(nmass)

# Model results
summary(nmass)
Anova(nmass)
r.squaredGLMM(nmass)

# Post-hoc tests
# Interaction between CO2 and n.trt
test(emtrends(nmass, pairwise~co2, "n.trt")) 
# Stronger stimulation in Nmass with increasing soil N under ambient CO2

# Interaction between inoc and n.trt
test(emtrends(nmass, pairwise~inoc, "n.trt"))
# Stronger stimulation in Nmass with increasing soil N in non-inoculated pots

# Individual effect of co2
emmeans(nmass, pairwise~co2)
# Nmass is larger in ambient CO2 treatment

# Individual effect of inoc
emmeans(nmass, pairwise~inoc)
# Nmass is larger in inoculated pots

# Individual effect of n.trt
test(emtrends(nmass, ~1, "n.trt"))
# Increasing N fertilization increases Nmass


##########################################################################
## Chlorophyll content
##########################################################################
chl.area <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
qqnorm(residuals(chl.area))
qqline(residuals(chl.area))
densityPlot(residuals(chl.area))
shapiro.test(residuals(chl.area))
outlierTest(chl.area)

# Model results
summary(chl.area)
Anova(chl.area)
r.squaredGLMM(chl.area)

# Post-hoc tests
# Interaction between CO2 and n.trt
test(emtrends(chl.area, pairwise~co2, "n.trt")) 
# Stronger stimulation in Chl.area with increasing soil N under ambient CO2

# Interaction between inoc and n.trt
test(emtrends(chl.area, pairwise~inoc, "n.trt"))
# Stronger stimulation in Chl.area with increasing soil N in non-inoculated pots

# Interaction between inoc and n.trt
emmeans(chl.area, pairwise~inoc*co2)
# Stronger stimulation in Chl.area with increasing soil N in non-inoculated pots


# Individual effect of co2
emmeans(chl.area, pairwise~co2)
# Chl.area is larger in ambient CO2 treatment

# Individual effect of inoc
emmeans(chl.area, pairwise~inoc)
# Chl.area is larger in inoculated pots

# Individual effect of n.trt
test(emtrends(chl.area, ~1, "n.trt"))
# Increasing N fertilization increases Chl.area


##########################################################################
## Vcmax25
##########################################################################
vcmax <- lmer(vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(vcmax)
qqnorm(residuals(vcmax))
qqline(residuals(vcmax))
densityPlot(residuals(vcmax))
shapiro.test(residuals(vcmax))
outlierTest(vcmax)

# Model results
summary(vcmax)
Anova(vcmax)
r.squaredGLMM(vcmax)

# Pairwise comparisons

# Two-way interaction between inoculation and n.trt
test(emtrends(vcmax, ~inoc, "n.trt"))
## Increasing N fertilization increased vcmax, but only in un-inoculated pots

# Individual effect of CO2
emmeans(vcmax, pairwise~co2)
## Elevated CO2 downregulated Vcmax

# Individual effect of inoculation
emmeans(vcmax, pairwise~inoc)
## Inoculated pots generally had larger Vcmax

# Individual effect of n.trt
test(emtrends(vcmax, ~1, "n.trt"))
## Vcmax generally increased with increasing N fertilization

##########################################################################
## Jmax25
##########################################################################
jmax <- lmer(jmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(jmax)
qqnorm(residuals(jmax))
qqline(residuals(jmax))
densityPlot(residuals(jmax))
shapiro.test(residuals(jmax))
outlierTest(jmax)

# Model results
summary(jmax)
Anova(jmax)
r.squaredGLMM(jmax)

# Pairwise comparisons

# Two-way interaction between inoculation and n.trt
test(emtrends(jmax, ~inoc, "n.trt"))
## Increasing N fertilization increased Jmax, but only in un-inoculated pots

# Individual effect of CO2
emmeans(jmax, pairwise~co2)
## Elevated CO2 downregulated Jmax

# Individual effect of inoculation
emmeans(jmax, pairwise~inoc)
## Inoculated pots generally had larger Jmax

# Individual effect of n.trt
test(emtrends(jmax, ~1, "n.trt"))
## Jmax generally increased with increasing N fertilization

##########################################################################
## Jmax25:Vcmax25
##########################################################################
df$jmax25.vcmax25[110] <- NA

jvmax <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(jvmax)
qqnorm(residuals(jvmax))
qqline(residuals(jvmax))
densityPlot(residuals(jvmax))
shapiro.test(residuals(jvmax))
outlierTest(jvmax)

# Model results
summary(jvmax)
Anova(jvmax)
r.squaredGLMM(jvmax)

# Pairwise comparisons

# Two-way interaction between inoculation and n.trt
test(emtrends(jvmax, ~inoc, "n.trt"))
## Increasing N fertilization increased Jmax:Vcmax, but only in 
## un-inoculated pots

# Individual effect of CO2
emmeans(jvmax, pairwise~co2)
## Elevated CO2 increases Jmax:Vcmax

# Individual effect of inoculation
emmeans(jvmax, pairwise~inoc)
## Non-inoculated pots generally had larger Jmax:Vcmax

# Individual effect of n.trt
test(emtrends(jvmax, ~1, "n.trt"))
## Jmax:Vcmax generally decreased with increasing N fertilization

##########################################################################
## Anet
##########################################################################
df$anet[110] <- NA

anet <- lmer(anet ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(anet)
qqnorm(residuals(anet))
qqline(residuals(anet))
densityPlot(residuals(anet))
shapiro.test(residuals(anet))
outlierTest(anet)

# Model results
summary(anet)
Anova(anet)
r.squaredGLMM(anet)

# Pairwise comparisons

# Two-way interaction between inoculation and n.trt
test(emtrends(anet, ~inoc, "n.trt"))
## Increasing N fertilization increased Anet in noninoculated pots; marginal
## negative effect of N fertilization in inoculated pots

# Individual effect of CO2
emmeans(anet, pairwise~co2)
## Elevated CO2 increases Anet

# Individual effect of inoculation
emmeans(anet, pairwise~inoc)
## Inoculated pots generally had larger Anet

# Individual effect of n.trt
test(emtrends(anet, ~1, "n.trt"))
## Anet generally decreased with increasing N fertilization


##########################################################################
## gsw
##########################################################################
df$gsw[78] <- NA

gsw <- lmer(gsw ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(gsw)
qqnorm(residuals(gsw))
qqline(residuals(gsw))
densityPlot(residuals(gsw))
shapiro.test(residuals(gsw))
outlierTest(gsw)

# Model results
summary(gsw)
Anova(gsw)
r.squaredGLMM(gsw)

# Pairwise comparisons

# Two-way interaction between inoculation and n.trt
test(emtrends(gsw, ~inoc, "n.trt"))
## Increasing N fertilization increased gs in noninoculated pots; decreased
## gsw with increasing fertilization in inoculated pots

# Marginal two-way interaction between co2 and n.trt
test(emtrends(gsw, ~co2, "n.trt"))
## Increasing N fertilization marginally increases gsw under ambient CO2;
## no N fertilization effect on gsw on elevated CO2

# Individual effect of CO2
emmeans(gsw, pairwise~co2)
## Elevated CO2 decreases stomatal conductance

# Individual effect of inoculation
emmeans(gsw, pairwise~inoc)
## Inoculated pots generally had higher stomatal conductance

# Individual effect of n.trt
test(emtrends(gsw, ~1, "n.trt"))
## N fertilization has no effect on stomatal conductance, although
## has positive trend

##########################################################################
## ci:ca
##########################################################################
cica <- lmer(ci.ca ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(cica)
qqnorm(residuals(cica))
qqline(residuals(cica))
densityPlot(residuals(cica))
shapiro.test(residuals(cica))
outlierTest(cica)

# Model results
summary(cica)
Anova(cica)
r.squaredGLMM(cica)

# Pairwise comparisons
## Only detectable pattern is a one-way effect of N fertilization on
## Ci:Ca:
test(emtrends(cica, ~1, "n.trt"))

##########################################################################
## stomatal limitation
##########################################################################
stomlim <- lmer(stomlim ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(stomlim)
qqnorm(residuals(stomlim))
qqline(residuals(stomlim))
densityPlot(residuals(stomlim))
shapiro.test(residuals(stomlim))
outlierTest(stomlim)

# Model results
summary(stomlim)
Anova(stomlim)
r.squaredGLMM(stomlim)

# Pairwise comparisons
## Only detectable pattern is a one-way effect of N fertilization on
## stomatal limitation:
test(emtrends(stomlim, ~1, "n.trt"))
## Increasing soil N increases stomatal limitation of photosynthesis

## and a marginal effect of inoculation on stomatal limitation:
emmeans(stomlim, pairwise~inoc)
## Inoculated pots have marginally higher stomatal limitation of photosynthesis
## than uninoculated pots

##########################################################################
## PNUE
##########################################################################
df$pnue[50] <- NA

pnue <- lmer(pnue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(pnue)
qqnorm(residuals(pnue))
qqline(residuals(pnue))
densityPlot(residuals(pnue))
shapiro.test(residuals(pnue))
outlierTest(pnue)

# Model results
summary(pnue)
Anova(pnue)
r.squaredGLMM(pnue)

# Pairwise comparisons

## Two-way interaction between inoculation and n.trt
test(emtrends(pnue, ~inoc, "n.trt"))
## Inoculated pots exhibited strong negative effect of increasing soil N
## on PNUE; negative but nonsignificant effect of increasing soil N on PNUE
## in non-inoculated pots

## Individual effect of CO2 on PNUE
emmeans(pnue, pairwise~co2)
## Ambient CO2 had lower PNUE than elevated CO2

## Individual effect of inoculation
emmeans(pnue, pairwise~inoc)
## Inoculated pots generally had greater PNUE than uninoculated pots

## Individual effect of soil N
test(emtrends(pnue, ~1, "n.trt"))
## Increasing soil N fertilization generally decreased PNUE

##########################################################################
## iWUE
##########################################################################
iwue <- lmer(iwue ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(iwue)
qqnorm(residuals(iwue))
qqline(residuals(iwue))
densityPlot(residuals(iwue))
shapiro.test(residuals(iwue))
outlierTest(iwue)

# Model results
summary(iwue)
Anova(iwue)
r.squaredGLMM(iwue)

# Pairwise comparisons

## Individual effect of n.trt on iWUE
test(emtrends(iwue, ~1, "n.trt"))
## Increasing soil N fertilization increases iWUE (interesting!!!!!)


##########################################################################
## Narea:gs
##########################################################################
narea.gs <- lmer(log(narea.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(narea.gs)
qqnorm(residuals(narea.gs))
qqline(residuals(narea.gs))
densityPlot(residuals(narea.gs))
shapiro.test(residuals(narea.gs))
outlierTest(narea.gs)

# Model results
summary(narea.gs)
Anova(narea.gs)
r.squaredGLMM(narea.gs)

# Pairwise comparisons
## Two-way interaction between inoculation and soil N
test(emtrends(narea.gs, ~inoc, "n.trt"))
## Increasing soil N increases Narea.gs in inoculated pots, but not uninoculated
## pots (although does trend in positive direction at p=0.146)

## Individual effect of CO2
emmeans(narea.gs, pairwise~co2)
## Narea.gs is greater under ambient CO2 

## Individual effect of inoculation
emmeans(narea.gs, pairwise~inoc)
## Narea.gs is greater in uninoculated pots

## Individual effect of n.trt on Narea:gs
test(emtrends(narea.gs, ~1, "n.trt"))
## Increasing soil N fertilization increases Narea:gs.

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- lmer(log(vcmax.gs) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(vcmax.gs)
qqnorm(residuals(vcmax.gs))
qqline(residuals(vcmax.gs))
densityPlot(residuals(vcmax.gs))
shapiro.test(residuals(vcmax.gs))
outlierTest(vcmax.gs)

# Model results
summary(vcmax.gs)
Anova(vcmax.gs)
r.squaredGLMM(vcmax.gs)

# Pairwise comparisons
## Individual effect of n.trt on Vcmax:gs
test(emtrends(vcmax.gs, ~1, "n.trt"))
## Increasing soil N fertilization increases Vcmax:gs.

##########################################################################
## Proportion of N in photosynthesis
##########################################################################
df$p.photo[50] <- NA
df$p.photo[df$p.photo > 1] <- NA

p.photo <- lmer(p.photo ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(p.photo)
qqnorm(residuals(p.photo))
qqline(residuals(p.photo))
densityPlot(residuals(p.photo))
shapiro.test(residuals(p.photo))
outlierTest(p.photo)

# Model results
summary(p.photo)
Anova(p.photo)
r.squaredGLMM(p.photo)

# Pairwise comparisons
## Two-way interaction between inoculation and n.trt
test(emtrends(p.photo, pairwise~inoc, "n.trt"))
## Negative effect of soil N on p.photo is stronger in inoculated pots
## than uninoculated pots. In fact, there is no effect of soil N on 
## photo in uninoculated pots (!?!), suggesting that deviations from 
## constant leaf N:photosynthesis relationships may only be apparent
## when there is no N limitation?????!?!?!?!? This seems important; added
## some detail to Notion notes page

## Individual effect of co2 on p.photo
emmeans(p.photo, pairwise~co2)
## Elevated CO2 generally has higher p.photo

## Individual effect of inoc on p.photo
emmeans(p.photo, pairwise~inoc)
## Inoculated pots generally have higher p.photo

## Individual effect of n.trt on p.photo
test(emtrends(p.photo, ~1, "n.trt"))
## Increasing soil N fertilization decreases p.photo

library(sjPlot)

p.pho.pred <- data.frame(get_model_data(p.photo, type = "pred", 
                                        terms = c("n.trt", "co2", "inoc")))


p.photo.inoc <- ggplot(data = subset(df, inoc == "inoc"), 
                       aes(x = n.trt, y = p.photo)) +
  geom_ribbon(data = subset(p.pho.pred, group_col == "amb" & facet == "inoc"), 
              aes(x = x, y = predicted, ymin = conf.low, 
                  ymax = conf.high), alpha = 0.25, fill = "#1965B0") +
  geom_line(data = subset(p.pho.pred, group_col == "amb" & facet == "inoc"), 
            aes(x = x, y = predicted), size = 1, color = "#1965B0") +
  geom_ribbon(data = subset(p.pho.pred, group_col == "elv" & facet == "inoc"), 
              aes(x = x, y = predicted, ymin = conf.low, 
                  ymax = conf.high), alpha = 0.25, fill = "#DC050C") +
  geom_line(data = subset(p.pho.pred, group_col == "elv" & facet == "inoc"), 
            aes(x = x, y = predicted), size = 1, color = "#DC050C") +
  geom_point(aes(fill = co2),
             size = 4, shape = 21, alpha = 0.75) +
  scale_fill_manual(values = c("#1965B0", "#DC050C"),
                    labels = c("Ambient", "Elevated")) +
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2)) +
  labs(title = "Inoculated",
       x = expression(bold("Soil N fertilization (ppm)")),
       y = expression(bold(rho["photo"]*"(g g"^"-1"*")")),
       fill = expression(bold("CO"["2"]*" treatment"))) +
  theme_bw(base_size = 18)

p.photo.noinoc <- ggplot(data = subset(df, inoc == "no.inoc"), 
                         aes(x = n.trt, y = p.photo)) +
  geom_ribbon(data = subset(p.pho.pred, group_col == "amb" & facet == "no.inoc"), 
              aes(x = x, y = predicted, ymin = conf.low, 
                  ymax = conf.high), alpha = 0.25, fill = "#1965B0") +
  geom_line(data = subset(p.pho.pred, group_col == "amb" & facet == "no.inoc"), 
            aes(x = x, y = predicted), size = 1, color = "#1965B0") +
  geom_ribbon(data = subset(p.pho.pred, group_col == "elv" & facet == "no.inoc"), 
              aes(x = x, y = predicted, ymin = conf.low, 
                  ymax = conf.high), alpha = 0.25, fill = "#DC050C") +
  geom_line(data = subset(p.pho.pred, group_col == "elv" & facet == "no.inoc"), 
            aes(x = x, y = predicted), size = 1, color = "#DC050C") +
  geom_point(aes(fill = co2),
             size = 4, shape = 23, alpha = 0.75) +
  scale_fill_manual(values = c("#1965B0", "#DC050C"),
                    labels = c("Ambient", "Elevated")) +
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2)) +
  labs(title = "Not inoculated",
       x = expression(bold("Soil N fertilization (ppm)")),
       y = expression(bold(rho["photo"]*"(g g"^"-1"*")")),
       fill = expression(bold("CO"["2"]*" treatment"))) +
  theme_bw(base_size = 18)


png("../working_drafts/figs/NxCO2xI_p.photo.fig.png",
    width = 10, height = 4, units = "in", res = 600)
ggarrange(p.photo.noinoc, p.photo.inoc, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()




##########################################################################
## Proportion of N in structure
##########################################################################
df$p.photo[50] <- NA

p.structure <- lmer(log(p.structure) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(p.structure)
qqnorm(residuals(p.structure))
qqline(residuals(p.structure))
densityPlot(residuals(p.structure))
shapiro.test(residuals(p.structure))
outlierTest(p.structure)

# Model results
summary(p.structure)
Anova(p.structure)
r.squaredGLMM(p.structure)

# Pairwise comparisons
## Two-way interaction between inoculation and n.trt
test(emtrends(p.structure, pairwise~inoc, "n.trt"))
## Negative effect of soil N on p.structure is greater in uninoculated 
## pots

## Individual effect of co2
emmeans(p.structure, pairwise~co2)
## Elevated CO2 generally has higher p.structure

## Individual effect of inoc
emmeans(p.structure, pairwise~inoc)
## Uninoculated pots generally have higher p.photo

## Individual effect of n.trt
test(emtrends(p.structure, ~1, "n.trt"))
## Increasing soil N fertilization decreases p.structure


##########################################################################
## Ncost
##########################################################################
ncost <- lmer(log(ncost) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(ncost)
qqnorm(residuals(ncost))
qqline(residuals(ncost))
densityPlot(residuals(ncost))
shapiro.test(residuals(ncost))
outlierTest(ncost)

# Model results
summary(ncost)
Anova(ncost)
r.squaredGLMM(ncost)

# Pairwise comparisons
## Two-way interaction between CO2 and soil N
test(emtrends(ncost, pairwise~co2, "n.trt"))
## Negative effect of increasing soil N is greater in ambient CO2

## Two-way interaction between inoc and soil N
test(emtrends(ncost, pairwise~inoc, "n.trt"))
## Negativ effect of increasing soil N is greater in non-inoculated pots

## Two-way interaction between CO2 and inoc
cld(emmeans(ncost, pairwise~co2*inoc))
## There is a stronger %difference between inoculation status at elevated
## CO2 than at ambient CO2. ## Specifically, noninoculated pots grown under 
## ambient CO2 had 9.5% higher Ncost, while noninoculated pots grown under
## elevated CO2 had 27.6% higher Ncost. Perhaps due to higher Ndemand?

## Individual effect of CO2
emmeans(ncost, pairwise~co2)
## elevated co2 generally has higher Ncost

## Individual effect of inoc
emmeans(ncost, pairwise~inoc)
## Noninoculated pots generally had higher Ncost

## Individual effect of soil N
test(emtrends(ncost, ~1, "n.trt"))
## Increasing soil N decreases Ncost

##########################################################################
## Total leaf area
##########################################################################
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(tla)
qqnorm(residuals(tla))
qqline(residuals(tla))
densityPlot(residuals(tla))
shapiro.test(residuals(tla))
outlierTest(tla)

# Model results
summary(tla)
Anova(tla)
r.squaredGLMM(tla)

# Pairwise comparisons
## Two-way interaction between CO2 and soil N
test(emtrends(tla, pairwise~co2, "n.trt"))
## Positive effect of increasing soil N is greater in elevated CO2

## Two-way interaction between inoc and soil N
test(emtrends(tla, pairwise~inoc, "n.trt"))
## Positive effect of increasing soil N is greater in non-inoculated pots

## Two-way interaction between CO2 and inoc
cld(emmeans(tla, pairwise~co2*inoc))
## There is a stronger %difference between inoculation status at elevated
## CO2 than at ambient CO2. ## Specifically, noninoculated pots grown under 
## ambient CO2 had 36.5% higher TLA, while noninoculated pots grown under
## elevated CO2 had 46.3% higher TLA. Perhaps due to higher Ndemand?

## Individual effect of CO2
emmeans(tla, pairwise~co2)
## elevated co2 generally has higher TLA

## Individual effect of inoc
emmeans(tla, pairwise~inoc)
## Inoculated pots generally had higher TLA

## Individual effect of soil N
test(emtrends(tla, ~1, "n.trt"))
## Increasing soil N decreases TLA

##########################################################################
## Total biomass
##########################################################################
tbio <- lmer(log(total.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(tbio)
qqnorm(residuals(tbio))
qqline(residuals(tbio))
densityPlot(residuals(tbio))
shapiro.test(residuals(tbio))
outlierTest(tbio)

# Model results
summary(tbio)
Anova(tbio)
r.squaredGLMM(tbio)

# Pairwise comparisons
## Two-way interaction between CO2 and soil N (marginal)
test(emtrends(tbio, pairwise~co2, "n.trt"))
## Marginally stronger positive effect of increasing soil N in elevated CO2

## Two-way interaction between inoc and soil N
test(emtrends(tbio, pairwise~inoc, "n.trt"))
## Positive effect of increasing soil N is greater in non-inoculated pots

## Individual effect of CO2
emmeans(tbio, pairwise~co2)
## elevated co2 generally has higher total biomass

## Individual effect of inoc
emmeans(tbio, pairwise~inoc)
## Inoculated pots generally had higher total biomass

## Individual effect of soil N
test(emtrends(tbio, ~1, "n.trt"))
## Increasing soil N decreases total biomass





