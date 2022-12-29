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
  mutate(n.trt = as.numeric(n.trt),
         rd25.vcmax25 = rd25 / vcmax25) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nodule.biomass < 0.05))
  ## filter all uninoculated pots that have nod biomass > 0.05 g

df.removed <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  dplyr::filter(inoc == "no.inoc" & nodule.biomass > 0.05)

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
test(emtrends(narea, pairwise~inoc, "n.trt"))
cld(emmeans(narea, pairwise~inoc*co2))

# Individual effects
emmeans(narea, pairwise~co2)
emmeans(narea, pairwise~inoc)
test(emtrends(narea, ~1, "n.trt"))

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
test(emtrends(marea, pairwise~co2, "n.trt")) 
test(emtrends(marea, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(marea, pairwise~co2)
emmeans(marea, pairwise~inoc)
test(emtrends(marea, ~1, "n.trt"))

##########################################################################
## Nmass
##########################################################################
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
test(emtrends(nmass, pairwise~co2, "n.trt")) 
test(emtrends(nmass, pairwise~inoc, "n.trt"))

# Individual effect of co2
emmeans(nmass, pairwise~co2)
emmeans(nmass, pairwise~inoc)
test(emtrends(nmass, ~1, "n.trt"))

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
test(emtrends(chl.area, pairwise~inoc, "n.trt"))
test(emtrends(chl.area, pairwise~co2, "n.trt")) 
test(emtrends(chl.area, pairwise~co2*inoc, "n.trt"))
emmeans(chl.area, pairwise~inoc*co2)

# Individual effect of co2
emmeans(chl.area, pairwise~co2)
emmeans(chl.area, pairwise~inoc)
test(emtrends(chl.area, ~1, "n.trt"))

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
test(emtrends(vcmax, ~inoc, "n.trt"))

# Individual effect of CO2
emmeans(vcmax, pairwise~co2)
emmeans(vcmax, pairwise~inoc)
test(emtrends(vcmax, ~1, "n.trt"))

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
test(emtrends(jmax, ~inoc, "n.trt"))

# Individual effects
emmeans(jmax, pairwise~co2)
emmeans(jmax, pairwise~inoc)
test(emtrends(jmax, ~1, "n.trt"))

##########################################################################
## Jmax25:Vcmax25
##########################################################################
df$jmax25.vcmax25[101] <- NA
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
test(emtrends(jvmax, pairwise~inoc, "n.trt"))
test(emtrends(jvmax, pairwise~co2, "n.trt"))

# Individual effects
emmeans(jvmax, pairwise~co2)
emmeans(jvmax, pairwise~inoc)
test(emtrends(jvmax, ~1, "n.trt"))

##########################################################################
## Anet
##########################################################################
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
test(emtrends(anet, ~inoc, "n.trt"))

# Individual effects
emmeans(anet, pairwise~co2)
emmeans(anet, pairwise~inoc)
test(emtrends(anet, ~1, "n.trt"))

##########################################################################
## Rd25
##########################################################################
df$rd25[df$rd25 < 0] <- NA

df$rd25[c(19, 29, 34, 57)] <- NA

rd25 <- lmer(rd25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(rd25)
qqnorm(residuals(rd25))
qqline(residuals(rd25))
densityPlot(residuals(rd25))
shapiro.test(residuals(rd25))
outlierTest(rd25)

# Model results
summary(rd25)
Anova(rd25)
r.squaredGLMM(rd25)

# Pairwise comparisons
test(emtrends(rd25, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(rd25, pairwise~co2)
emmeans(rd25, pairwise~inoc)
test(emtrends(rd25, ~1, "n.trt"))


##########################################################################
## gsw
##########################################################################
df$gsw[70] <- NA
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
test(emtrends(gsw, pairwise~inoc, "n.trt"))

# Individual effect
emmeans(gsw, pairwise~co2)
emmeans(gsw, pairwise~inoc)
test(emtrends(gsw, ~1, "n.trt"))

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
emmeans(stomlim, pairwise~co2*inoc)

# Individual effects
test(emtrends(stomlim, ~1, "n.trt"))
emmeans(stomlim, pairwise~inoc)

##########################################################################
## Proportion of N in photosynthesis
##########################################################################
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
test(emtrends(p.photo, pairwise~inoc, "n.trt"))
test(emtrends(p.photo, pairwise~co2*inoc, "n.trt"))


## Individual effects
emmeans(p.photo, pairwise~co2)
emmeans(p.photo, pairwise~inoc)
test(emtrends(p.photo, ~1, "n.trt"))

##########################################################################
## Proportion of N in structure
##########################################################################
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
test(emtrends(p.structure, pairwise~co2, "n.trt"))
test(emtrends(p.structure, pairwise~inoc, "n.trt"))
emmeans(p.structure, pairwise~co2*inoc, type = "response")

## Individual effects
emmeans(p.structure, pairwise~co2, type = "response")
emmeans(p.structure, pairwise~inoc, type = "response")
test(emtrends(p.structure, ~1, "n.trt"))

##########################################################################
## Proportion of N in Rubisco
##########################################################################
df$p.rubisco[df$p.rubisco>1] <- NA
p.rub <- lmer(p.rubisco ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(p.rub)
qqnorm(residuals(p.rub))
qqline(residuals(p.rub))
densityPlot(residuals(p.rub))
shapiro.test(residuals(p.rub))
outlierTest(p.rub)

# Model results
summary(p.rub)
Anova(p.rub)
r.squaredGLMM(p.rub)

# Pairwise comparisons
test(emtrends(p.rub, pairwise~inoc, "n.trt"))

## Individual effect
emmeans(p.rub, pairwise~co2)
emmeans(p.rub, pairwise~inoc)
test(emtrends(p.rub, ~1, "n.trt"))

##########################################################################
## Proportion of N in bioenergetics
##########################################################################
df$p.bioe[c(45)] <- NA
p.bioe <- lmer(p.bioe ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(p.bioe)
qqnorm(residuals(p.bioe))
qqline(residuals(p.bioe))
densityPlot(residuals(p.bioe))
shapiro.test(residuals(p.bioe))
outlierTest(p.bioe)

# Model results
summary(p.bioe)
Anova(p.bioe)
r.squaredGLMM(p.bioe)

# Pairwise comparisons
test(emtrends(p.bioe, pairwise~inoc, "n.trt"))

## Individual effects
emmeans(p.bioe, pairwise~co2)
emmeans(p.bioe, pairwise~inoc)
test(emtrends(p.bioe, ~1, "n.trt"))

##########################################################################
## Proportion of N in light harvesting
##########################################################################
df$p.lightharv[c(25, 39, 45)] <- NA
p.light <- lmer(p.lightharv ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(p.light)
qqnorm(residuals(p.light))
qqline(residuals(p.light))
densityPlot(residuals(p.light))
shapiro.test(residuals(p.light))
outlierTest(p.light)

# Model results
summary(p.light)
Anova(p.light)
r.squaredGLMM(p.light)

# Pairwise comparisons
test(emtrends(p.light, pairwise~co2*inoc, "n.trt"))
test(emtrends(p.light, pairwise~inoc, "n.trt"))

## Individual effects
emmeans(p.light, pairwise~co2)
emmeans(p.light, pairwise~inoc)
test(emtrends(p.light, ~1, "n.trt"))

##########################################################################
## PNUE
##########################################################################
df$pnue[45] <- NA
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
test(emtrends(pnue, ~inoc, "n.trt"))


## Individual effects
emmeans(pnue, pairwise~co2)
emmeans(pnue, pairwise~inoc)
test(emtrends(pnue, ~1, "n.trt"))

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
test(emtrends(narea.gs, ~inoc, "n.trt"))

## Individual effects
emmeans(narea.gs, pairwise~co2)
emmeans(narea.gs, pairwise~inoc)
test(emtrends(narea.gs, ~1, "n.trt"))

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
emmeans(vcmax.gs, pairwise~co2*inoc, 
        type = "response")

test(emtrends(vcmax.gs, ~1, "n.trt"))


##########################################################################
## Ncost
##########################################################################
df$ncost[c(101, 102)] <- NA

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
## Three-way interaction
test(emtrends(ncost, pairwise~co2*inoc, "n.trt"))

## Two-way interaction between CO2 and soil N
test(emtrends(ncost, pairwise~co2, "n.trt"))
test(emtrends(ncost, pairwise~inoc, "n.trt"))
emmeans(ncost, pairwise~co2*inoc)

## Individual effects
emmeans(ncost, pairwise~co2)
emmeans(ncost, pairwise~inoc)
test(emtrends(ncost, ~1, "n.trt"))

##########################################################################
## Belowground carbon biomass
##########################################################################
cbg <- lmer(log(cbg) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(cbg)
qqnorm(residuals(cbg))
qqline(residuals(cbg))
densityPlot(residuals(cbg))
shapiro.test(residuals(cbg))
outlierTest(cbg)

# Model results
summary(cbg)
Anova(cbg)
r.squaredGLMM(cbg)

# Pairwise comparisons
test(emtrends(cbg, pairwise~inoc, "n.trt"))
test(emtrends(cbg, pairwise~co2, "n.trt"))
emmeans(cbg, pairwise~inoc*co2, type = "response")

## Individual effects
emmeans(cbg, pairwise~co2, type = "response")
emmeans(cbg, pairwise~inoc, type = "response")
test(emtrends(cbg, ~1, "n.trt"))

##########################################################################
## Whole plant nitrogen
##########################################################################
wpn <- lmer(log(wpn) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(wpn)
qqnorm(residuals(wpn))
qqline(residuals(wpn))
densityPlot(residuals(wpn))
shapiro.test(residuals(wpn))
outlierTest(wpn)

# Model results
summary(wpn)
Anova(wpn)
r.squaredGLMM(wpn)

# Pairwise comparisons
test(emtrends(wpn, pairwise~inoc, "n.trt"))

## Individual effects
emmeans(wpn, pairwise~co2, type = "response")
emmeans(wpn, pairwise~inoc, type = "response")
test(emtrends(wpn, ~1, "n.trt", regrid = "response"))

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
test(emtrends(tla, pairwise~co2, "n.trt"))
test(emtrends(tla, pairwise~inoc, "n.trt"))
emmeans(tla, pairwise~co2*inoc)

## Individual effects
emmeans(tla, pairwise~co2)
emmeans(tla, pairwise~inoc)
test(emtrends(tla, ~1, "n.trt"))

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
test(emtrends(tbio, pairwise~inoc, "n.trt"))
test(emtrends(tbio, pairwise~co2, "n.trt"))
emmeans(tbio, pairwise~co2*inoc)

## Individual effects
emmeans(tbio, pairwise~co2)
emmeans(tbio, pairwise~inoc)
test(emtrends(tbio, ~1, "n.trt"))

##########################################################################
## Root nodule biomass: root biomass
##########################################################################
nod.bio <- lmer(sqrt(nodule.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(nod.bio)
qqnorm(residuals(nod.bio))
qqline(residuals(nod.bio))
densityPlot(residuals(nod.bio))
shapiro.test(residuals(nod.bio))
outlierTest(nod.bio)

# Model results
summary(nod.bio)
Anova(nod.bio)
r.squaredGLMM(nod.bio)

# Pairwise comparisons
test(emtrends(nod.bio, pairwise~inoc, "n.trt"))
test(emtrends(nod.bio, pairwise~co2, "n.trt"))

## Individual effects
emmeans(nod.bio, pairwise~co2)
emmeans(nod.bio, pairwise~inoc)
test(emtrends(nod.bio, ~1, "n.trt"))

##########################################################################
## Root nodule biomass: root biomass
##########################################################################
df$nod.root.ratio <- df$nodule.biomass / df$root.biomass
df$nod.root.ratio[df$nod.root.ratio > 0.05 & df$inoc == "no.inoc"] <- NA



nod.root.ratio <- lmer(sqrt(nod.root.ratio) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model assumptions
plot(nod.root.ratio)
qqnorm(residuals(nod.root.ratio))
qqline(residuals(nod.root.ratio))
densityPlot(residuals(nod.root.ratio))
shapiro.test(residuals(nod.root.ratio))
outlierTest(nod.root.ratio)

# Model results
summary(nod.root.ratio)
Anova(nod.root.ratio)
r.squaredGLMM(nod.root.ratio)

# Pairwise comparisons
## Three-way interaction between inoculation and N fertilization
test(emtrends(nod.root.ratio, pairwise~inoc*co2, "n.trt"))
test(emtrends(nod.root.ratio, pairwise~inoc, "n.trt"))
test(emtrends(nod.root.ratio, pairwise~co2, "n.trt"))

emmeans(nod.root.ratio, pairwise~co2*inoc)

## Individual effects
emmeans(nod.root.ratio, pairwise~co2)
emmeans(nod.root.ratio, pairwise~inoc)
test(emtrends(nod.root.ratio, ~1, "n.trt"))

##########################################################################
## %Ndfa
##########################################################################
# ndfa <- lmer(ndfa ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
# 
# # Check model assumptions
# plot(ndfa)
# qqnorm(residuals(ndfa))
# qqline(residuals(ndfa))
# densityPlot(residuals(ndfa))
# shapiro.test(residuals(ndfa))
# outlierTest(ndfa)
# 
# # Model results
# summary(ndfa)
# Anova(ndfa)
# r.squaredGLMM(ndfa)

# Pairwise comparisons

##########################################################################
## Table 1: Leaf N content
##########################################################################
narea.coefs <- data.frame(summary(narea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.narea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.narea) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.narea) %>%
  print(., row.names = FALSE)

narea.table <- data.frame(Anova(narea)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(narea.coefs) %>%
  dplyr::select(treatment, df = Df, coef.narea, 
                chisq.narea = Chisq, pval.narea = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.narea:pval.narea, round, digits = 3),
         chisq.narea = ifelse(chisq.narea <0.001 & chisq.narea >= 0, 
                              "<0.001", chisq.narea),
         pval.narea = ifelse(pval.narea <0.001 & pval.narea >= 0, 
                             "<0.001", pval.narea)) %>%
  arrange(treatment)

nmass.coefs <- data.frame(summary(nmass)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nmass = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.nmass) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.nmass) %>%
  print(., row.names = FALSE)

nmass.table <- data.frame(Anova(nmass)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nmass.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nmass, 
                chisq.nmass = Chisq, pval.nmass = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.nmass:pval.nmass, round, digits = 3),
         chisq.nmass = ifelse(chisq.nmass <0.001 & chisq.nmass >= 0, 
                              "<0.001", chisq.nmass),
         pval.nmass = ifelse(pval.nmass <0.001 & pval.nmass >= 0, 
                             "<0.001", pval.nmass)) %>%
  arrange(treatment)

marea.coefs <- data.frame(summary(marea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.marea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.marea) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.marea) %>%
  print(., row.names = FALSE)

marea.table <- data.frame(Anova(marea)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(marea.coefs) %>%
  dplyr::select(treatment, df = Df, coef.marea, 
                chisq.marea = Chisq, pval.marea = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.marea:pval.marea, round, digits = 3),
         chisq.marea = ifelse(chisq.marea <0.001 & chisq.marea >= 0, 
                              "<0.001", chisq.marea),
         pval.marea = ifelse(pval.marea <0.001 & pval.marea >= 0, 
                             "<0.001", pval.marea)) %>%
  arrange(treatment)

chl.area.coefs <- data.frame(summary(chl.area)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.chl.area = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.chl.area) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.chl.area) %>%
  print(., row.names = FALSE)

chl.area.table <- data.frame(Anova(chl.area)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(chl.area.coefs) %>%
  dplyr::select(treatment, df = Df, coef.chl.area, 
                chisq.chl.area = Chisq, pval.chl.area = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.chl.area:pval.chl.area, round, digits = 3),
         chisq.chl.area = ifelse(chisq.chl.area <0.001 & chisq.chl.area >= 0, 
                              "<0.001", chisq.chl.area),
         pval.chl.area = ifelse(pval.chl.area <0.001 & pval.chl.area >= 0, 
                             "<0.001", pval.chl.area)) %>%
  arrange(treatment)

table3 <- narea.table %>% full_join(nmass.table) %>% 
  full_join(marea.table) %>% full_join(chl.area.table) %>%
  replace(is.na(.), "-")
write.csv(table3, file = "../working_drafts/tables/NxCO2xI_table1_leafN.csv",
          row.names = FALSE)

##########################################################################
## Table 2: Gas exchange
##########################################################################
vcmax.coefs <- data.frame(summary(vcmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.vcmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.vcmax) %>%
  print(., row.names = FALSE)

vcmax.table <- data.frame(Anova(vcmax)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(vcmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.vcmax, 
                chisq.vcmax = Chisq, pval.vcmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.vcmax:pval.vcmax, round, digits = 3),
         chisq.vcmax = ifelse(chisq.vcmax < 0.001 & chisq.vcmax >= 0, 
                              "<0.001", chisq.vcmax),
         pval.vcmax = ifelse(pval.vcmax < 0.001 & pval.vcmax >= 0, 
                             "<0.001", pval.vcmax)) %>%
  arrange(treatment)

jmax.coefs <- data.frame(summary(jmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.jmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.jmax) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.jmax) %>%
  print(., row.names = FALSE)

jmax.table <- data.frame(Anova(jmax)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(jmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.jmax, 
                chisq.jmax = Chisq, pval.jmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.jmax:pval.jmax, round, digits = 3),
         chisq.jmax = ifelse(chisq.jmax < 0.001 & chisq.jmax >= 0, 
                              "<0.001", chisq.jmax),
         pval.jmax = ifelse(pval.jmax < 0.001 & pval.jmax >= 0, 
                             "<0.001", pval.jmax)) %>%
  arrange(treatment)

jvmax.coefs <- data.frame(summary(jvmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.jvmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.jvmax) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.jvmax) %>%
  print(., row.names = FALSE)

jvmax.table <- data.frame(Anova(jvmax)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(jvmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.jvmax, 
                chisq.jvmax = Chisq, pval.jvmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.jvmax:pval.jvmax, round, digits = 3),
         chisq.jvmax = ifelse(chisq.jvmax < 0.001 & chisq.jvmax >= 0, 
                             "<0.001", chisq.jvmax),
         pval.jvmax = ifelse(pval.jvmax < 0.001 & pval.jvmax >= 0, 
                            "<0.001", pval.jvmax)) %>%
  arrange(treatment)


rd25.coefs <- data.frame(summary(rd25)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.rd25 = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.rd25) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.rd25) %>%
  print(., row.names = FALSE)

rd25.table <- data.frame(Anova(rd25)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(rd25.coefs) %>%
  dplyr::select(treatment, df = Df, coef.rd25, 
                chisq.rd25 = Chisq, pval.rd25 = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.rd25:pval.rd25, round, digits = 3),
         chisq.rd25 = ifelse(chisq.rd25 < 0.001 & chisq.rd25 >= 0, 
                            "<0.001", chisq.rd25),
         pval.rd25 = ifelse(pval.rd25 < 0.001 & pval.rd25 >= 0, 
                           "<0.001", pval.rd25)) %>%
  arrange(treatment)

gsw.coefs <- data.frame(summary(gsw)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.gsw = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.gsw) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.gsw) %>%
  print(., row.names = FALSE)

gsw.table <- data.frame(Anova(gsw)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(gsw.coefs) %>%
  dplyr::select(treatment, df = Df, coef.gsw, 
                chisq.gsw = Chisq, pval.gsw = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.gsw:pval.gsw, round, digits = 3),
         chisq.gsw = ifelse(chisq.gsw < 0.001 & chisq.gsw >= 0, 
                              "<0.001", chisq.gsw),
         pval.gsw = ifelse(pval.gsw < 0.001 & pval.gsw >= 0, 
                             "<0.001", pval.gsw)) %>%
  arrange(treatment)

l.coefs <- data.frame(summary(stomlim)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.l = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.l) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.l) %>%
  print(., row.names = FALSE)

l.table <- data.frame(Anova(stomlim)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(l.coefs) %>%
  dplyr::select(treatment, df = Df, coef.l, 
                chisq.l = Chisq, pval.l = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.l:pval.l, round, digits = 3),
         chisq.l = ifelse(chisq.l < 0.001 & chisq.l >= 0, 
                             "<0.001", chisq.l),
         pval.l = ifelse(pval.l < 0.001 & pval.l >= 0, 
                            "<0.001", pval.l)) %>%
  arrange(treatment)

table4 <- vcmax.table %>% full_join(jmax.table) %>% full_join(jvmax.table) %>%
  full_join(rd25.table) %>% full_join(gsw.table) %>% full_join(l.table) %>%
  replace(is.na(.), "-")
write.csv(table4, file = "../working_drafts/tables/NxCO2xI_table2_gasEx.csv",
          row.names = FALSE)

##########################################################################
## Table 3: Prop leaf N to photosynthesis, structure, etc.
##########################################################################
p.rub.coefs <- data.frame(summary(p.rub)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.p.rub = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.rub) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.p.rub) %>%
  print(., row.names = FALSE)

p.rub.table <- data.frame(Anova(p.rub)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(p.rub.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.rub, 
                chisq.p.rub = Chisq, pval.p.rub = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.p.rub:pval.p.rub, round, digits = 3),
         chisq.p.rub = ifelse(chisq.p.rub < 0.001 & chisq.p.rub >= 0, 
                             "<0.001", chisq.p.rub),
         pval.p.rub = ifelse(pval.p.rub < 0.001 & pval.p.rub >= 0, 
                            "<0.001", pval.p.rub)) %>%
  arrange(treatment)

p.bioe.coefs <- data.frame(summary(p.bioe)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.p.bioe = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.bioe) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.p.bioe) %>%
  print(., row.names = FALSE)

p.bioe.table <- data.frame(Anova(p.bioe)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(p.bioe.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.bioe, 
                chisq.p.bioe = Chisq, pval.p.bioe = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.p.bioe:pval.p.bioe, round, digits = 3),
         chisq.p.bioe = ifelse(chisq.p.bioe < 0.001 & chisq.p.bioe >= 0, 
                              "<0.001", chisq.p.bioe),
         pval.p.bioe = ifelse(pval.p.bioe < 0.001 & pval.p.bioe >= 0, 
                             "<0.001", pval.p.bioe)) %>%
  arrange(treatment)

p.light.coefs <- data.frame(summary(p.light)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.p.light = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.light) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.p.light) %>%
  print(., row.names = FALSE)

p.light.table <- data.frame(Anova(p.light)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(p.light.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.light, 
                chisq.p.light = Chisq, pval.p.light = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.p.light:pval.p.light, round, digits = 3),
         chisq.p.light = ifelse(chisq.p.light < 0.001 & chisq.p.light >= 0, 
                               "<0.001", chisq.p.light),
         pval.p.light = ifelse(pval.p.light < 0.001 & pval.p.light >= 0, 
                              "<0.001", pval.p.light)) %>%
  arrange(treatment)

p.photo.coefs <- data.frame(summary(p.photo)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.p.photo = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.photo) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.p.photo) %>%
  print(., row.names = FALSE)

p.photo.table <- data.frame(Anova(p.photo)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(p.photo.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.photo, 
                chisq.p.photo = Chisq, pval.p.photo = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.p.photo:pval.p.photo, round, digits = 3),
         chisq.p.photo = ifelse(chisq.p.photo < 0.001 & chisq.p.photo >= 0, 
                               "<0.001", chisq.p.photo),
         pval.p.photo = ifelse(pval.p.photo < 0.001 & pval.p.photo >= 0, 
                              "<0.001", pval.p.photo)) %>%
  arrange(treatment)

p.str.coefs <- data.frame(summary(p.structure)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.p.str = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.p.str) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.p.str) %>%
  print(., row.names = FALSE)

p.str.table <- data.frame(Anova(p.structure)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(p.str.coefs) %>%
  dplyr::select(treatment, df = Df, coef.p.str, 
                chisq.p.str = Chisq, pval.p.str = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.p.str:pval.p.str, round, digits = 3),
         chisq.p.str = ifelse(chisq.p.str < 0.001 & chisq.p.str >= 0, 
                                "<0.001", chisq.p.str),
         pval.p.str = ifelse(pval.p.str < 0.001 & pval.p.str >= 0, 
                               "<0.001", pval.p.str)) %>%
  arrange(treatment)

table5 <- p.rub.table %>% full_join(p.bioe.table) %>% 
  full_join(p.light.table) %>% full_join(p.photo.table) %>%
  full_join(p.str.table) %>%
  replace(is.na(.), "-")

write.csv(table5, file = "../working_drafts/tables/NxCO2xI_table3_propN.csv",
          row.names = FALSE)

##########################################################################
## Table 4: PNUE/iWUE
##########################################################################
pnue.coefs <- data.frame(summary(pnue)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.pnue = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.pnue) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

pnue.table <- data.frame(Anova(pnue)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(pnue.coefs) %>%
  dplyr::select(treatment, df = Df, coef.pnue, 
                chisq.pnue = Chisq, pval.pnue = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.pnue:pval.pnue, round, digits = 3),
         chisq.pnue = ifelse(chisq.pnue < 0.001 & chisq.pnue >= 0, 
                              "<0.001", chisq.pnue),
         pval.pnue = ifelse(pval.pnue < 0.001 & pval.pnue >= 0, 
                             "<0.001", pval.pnue)) %>%
  arrange(treatment)

iwue.coefs <- data.frame(summary(iwue)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.iwue = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.iwue) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

iwue.table <- data.frame(Anova(iwue)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(iwue.coefs) %>%
  dplyr::select(treatment, df = Df, coef.iwue, 
                chisq.iwue = Chisq, pval.iwue = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.iwue:pval.iwue, round, digits = 3),
         chisq.iwue = ifelse(chisq.iwue < 0.001 & chisq.iwue >= 0, 
                             "<0.001", chisq.iwue),
         pval.iwue = ifelse(pval.iwue < 0.001 & pval.iwue >= 0, 
                            "<0.001", pval.iwue)) %>%
  arrange(treatment)

narea.gs.coefs <- data.frame(summary(narea.gs)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.narea.gs = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.narea.gs) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

narea.gs.table <- data.frame(Anova(narea.gs)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(narea.gs.coefs) %>%
  dplyr::select(treatment, df = Df, coef.narea.gs, 
                chisq.narea.gs = Chisq, pval.narea.gs = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.narea.gs:pval.narea.gs, round, digits = 3),
         chisq.narea.gs = ifelse(chisq.narea.gs < 0.001 & chisq.narea.gs >= 0, 
                             "<0.001", chisq.narea.gs),
         pval.narea.gs = ifelse(pval.narea.gs < 0.001 & pval.narea.gs >= 0, 
                            "<0.001", pval.narea.gs)) %>%
  arrange(treatment)

vcmax.gs.coefs <- data.frame(summary(vcmax.gs)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.vcmax.gs = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax.gs) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

vcmax.gs.table <- data.frame(Anova(vcmax.gs)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(vcmax.gs.coefs) %>%
  dplyr::select(treatment, df = Df, coef.vcmax.gs, 
                chisq.vcmax.gs = Chisq, pval.vcmax.gs = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.vcmax.gs:pval.vcmax.gs, round, digits = 3),
         chisq.vcmax.gs = ifelse(chisq.vcmax.gs < 0.001 & chisq.vcmax.gs >= 0, 
                                 "<0.001", chisq.vcmax.gs),
         pval.vcmax.gs = ifelse(pval.vcmax.gs < 0.001 & pval.vcmax.gs >= 0, 
                                "<0.001", pval.vcmax.gs)) %>%
  arrange(treatment)

table6 <- pnue.table %>% full_join(iwue.table) %>% full_join(narea.gs.table) %>%
  full_join(vcmax.gs.table) %>%
  replace(is.na(.), "-")
write.csv(table6, file = "../working_drafts/tables/NxCO2xI_table4_pnue_iwue.csv",
          row.names = FALSE)


##########################################################################
## Table 5: Whole plant traits
##########################################################################
ncost.coefs <- data.frame(summary(ncost)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.ncost = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.ncost) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

ncost.table <- data.frame(Anova(ncost)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(ncost.coefs) %>%
  dplyr::select(treatment, df = Df, coef.ncost, 
                chisq.ncost = Chisq, pval.ncost = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.ncost:pval.ncost, round, digits = 3),
         chisq.ncost = ifelse(chisq.ncost <0.001 & chisq.ncost >= 0, 
                              "<0.001", chisq.ncost),
         pval.ncost = ifelse(pval.ncost <0.001 & pval.ncost >= 0, 
                             "<0.001", pval.ncost)) %>%
  arrange(treatment)

cbg.coefs <- data.frame(summary(cbg)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.cbg = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.cbg) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

cbg.table <- data.frame(Anova(cbg)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(cbg.coefs) %>%
  dplyr::select(treatment, df = Df, coef.cbg, 
                chisq.cbg = Chisq, pval.cbg = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.cbg:pval.cbg, round, digits = 3),
         chisq.cbg = ifelse(chisq.cbg <0.001 & chisq.cbg >= 0, 
                            "<0.001", chisq.cbg),
         pval.cbg = ifelse(pval.cbg <0.001 & pval.cbg >= 0, 
                           "<0.001", pval.cbg)) %>%
  arrange(treatment)

wpn.coefs <- data.frame(summary(wpn)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.wpn = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.wpn) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

wpn.table <- data.frame(Anova(wpn)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(wpn.coefs) %>%
  dplyr::select(treatment, df = Df, coef.wpn, 
                chisq.wpn = Chisq, pval.wpn = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.wpn:pval.wpn, round, digits = 3),
         chisq.wpn = ifelse(chisq.wpn <0.001 & chisq.wpn >= 0, 
                            "<0.001", chisq.wpn),
         pval.wpn = ifelse(pval.wpn <0.001 & pval.wpn >= 0, 
                           "<0.001", pval.wpn)) %>%
  arrange(treatment)

tla.coefs <- data.frame(summary(tla)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.tla = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.tla) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

tla.table <- data.frame(Anova(tla)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(tla.coefs) %>%
  dplyr::select(treatment, df = Df, coef.tla, 
                chisq.tla = Chisq, pval.tla = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.tla:pval.tla, round, digits = 3),
         chisq.tla = ifelse(chisq.tla <0.001 & chisq.tla >= 0, 
                            "<0.001", chisq.tla),
         pval.tla = ifelse(pval.tla <0.001 & pval.tla >= 0, 
                           "<0.001", pval.tla)) %>%
  arrange(treatment)

tbio.coefs <- data.frame(summary(tbio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.tbio = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.tbio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

tbio.table <- data.frame(Anova(tbio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(tbio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.tbio, 
                chisq.tbio = Chisq, pval.tbio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.tbio:pval.tbio, round, digits = 3),
         chisq.tbio = ifelse(chisq.tbio < 0.001 & chisq.tbio >= 0, 
                             "<0.001", chisq.tbio),
         pval.tbio = ifelse(pval.tbio < 0.001 & pval.tbio >= 0, 
                            "<0.001", pval.tbio)) %>%
  arrange(treatment)



table5 <- ncost.table %>% full_join(cbg.table) %>% full_join(wpn.table) %>%
  full_join(tla.table) %>% full_join(tbio.table) %>%
  replace(is.na(.), "-")
write.csv(table5, file = "../working_drafts/tables/NxCO2xI_table5_WP.csv",
          row.names = FALSE)

##########################################################################
## Table 6: Nitrogen fixation
##########################################################################
nodroot.coefs <- data.frame(summary(nod.root.ratio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nodroot = format(Estimate, scientific = TRUE, digits = 3),
         se.nodroot = round(Std..Error, digits = 3),
         t.value.nodroot = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.nodroot, se.nodroot, t.value.nodroot) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.nodroot) %>%
  print(., row.names = FALSE)

nodroot.table <- data.frame(Anova(nod.root.ratio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nodroot.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nodroot, 
                chisq.nodroot = Chisq, pval.nodroot = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.nodroot:pval.nodroot, round, digits = 3),
         chisq.nodroot = ifelse(chisq.nodroot <0.001 & chisq.nodroot >= 0, 
                                "<0.001", chisq.nodroot),
         pval.nodroot = ifelse(pval.nodroot <0.001 & pval.nodroot >= 0, 
                               "<0.001", pval.nodroot)) %>%
  arrange(treatment)

nodbio.coefs <- data.frame(summary(nod.bio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nodbio = format(Estimate, scientific = TRUE, digits = 3),
         se.nodbio = round(Std..Error, digits = 3),
         t.value.nodbio = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.nodbio, se.nodbio, t.value.nodbio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.nodbio) %>%
  print(., row.names = FALSE)

nodbio.table <- data.frame(Anova(nod.bio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nodbio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nodbio, 
                chisq.nodbio = Chisq, pval.nodbio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.nodbio:pval.nodbio, round, digits = 3),
         chisq.nodbio = ifelse(chisq.nodbio <0.001 & chisq.nodbio >= 0, 
                               "<0.001", chisq.nodbio),
         pval.nodbio = ifelse(pval.nodbio <0.001 & pval.nodbio >= 0, 
                              "<0.001", pval.nodbio)) %>%
  arrange(treatment)

# ndfa.coefs <- data.frame(summary(ndfa)$coefficient) %>%
#   mutate(treatment = row.names(.),
#          coef.ndfa = format(Estimate, scientific = TRUE, digits = 3),
#          se.ndfa = round(Std..Error, digits = 3),
#          t.value.ndfa = round(t.value, digits = 3)) %>%
#   dplyr::select(treatment, coef.ndfa, se.ndfa, t.value.ndfa) %>%
#   mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
#                        "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
#   dplyr::select(treatment, coef.ndfa) %>%
#   print(., row.names = FALSE)
# 
# ndfa.table <- data.frame(Anova(ndfa)) %>%
#   mutate(treatment = row.names(.)) %>%
#   full_join(ndfa.coefs) %>%
#   dplyr::select(treatment, df = Df, coef.ndfa, 
#                 chisq.ndfa = Chisq, pval.ndfa = Pr..Chisq.) %>%
#   mutate(treatment = factor(treatment, 
#                             levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
#                                        "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
#          across(chisq.ndfa:pval.ndfa, round, digits = 3),
#          chisq.ndfa = ifelse(chisq.ndfa < 0.001 & chisq.ndfa >= 0, 
#                                "<0.001", chisq.ndfa),
#          pval.ndfa = ifelse(pval.ndfa < 0.001 & pval.ndfa >= 0, 
#                               "<0.001", pval.ndfa)) %>%
#   arrange(treatment)

table6 <- nodbio.table %>%
   full_join(nodroot.table) %>% 
  # full_join(ndfa.table) %>%
  replace(is.na(.), "-")
write.csv(table6, file = "../working_drafts/tables/NxCO2xI_table6_nfix.csv",
          row.names = FALSE)



