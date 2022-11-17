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
  mutate(n.trt = as.numeric(n.trt)) %>%
  dplyr::filter(id != "a_n_0_109" & id != "a_n_0_112" & id != "e_n_70_47" &
           id != "a_n_70_118" & id != "a_n_105_122" & id != "a_n_280_135")

##########################################################################
## Ncost
##########################################################################
df$ncost[c(108, 109)] <- NA

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
test(emtrends(ncost, ~1, "n.trt", regrid = "response"))

## Two-way interaction between CO2 and soil N
cld(emtrends(ncost, pairwise~co2*inoc, "n.trt"))
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
## Belowground carbon biomass
##########################################################################
cbg <- lmer(cbg ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
## Two-way interaction between CO2 and soil N
test(emtrends(cbg, pairwise~inoc, "n.trt"))
## Positive effect of increasing soil N is greater in uninoculated pots

## Two-way interaction between inoc and soil N
test(emtrends(cbg, pairwise~co2, "n.trt"))
## Negativ effect of increasing soil N is marginally greater in elevated CO2

## Two-way interaction between CO2 and inoc
cld(emmeans(cbg, pairwise~inoc*co2, type = "response"))
## There is a stronger %difference between inoculation status at ambient
## CO2 than at elevated CO2.

## Individual effect of CO2
emmeans(cbg, pairwise~co2)
## elevated co2 generally has higher values

## Individual effect of inoc
emmeans(cbg, pairwise~inoc)
## Inoculated pots generally had higher values

## Individual effect of soil N
test(emtrends(cbg, ~1, "n.trt", regrid = "response"))
emmeans(cbg, ~1, "n.trt", at = list(n.trt = 0), type = "response")
## Increasing soil N increases cbg
## Eq: 0.0023x + 0.469 

##########################################################################
## Whole plant nitrogen
##########################################################################
wpn <- lmer(sqrt(wpn) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
## Two-way interaction between CO2 and soil N
test(emtrends(wpn, pairwise~inoc, "n.trt"))
## Positive effect of increasing soil N is greater in uninoculated pots

## Two-way interaction between inoc and soil N
test(emtrends(wpn, pairwise~co2, "n.trt"))
## Negativ effect of increasing soil N is greater in non-inoculated pots

## Two-way interaction between CO2 and inoc
cld(emmeans(wpn, pairwise~inoc*co2, type = "response"))
## There is a stronger %difference between inoculation status at elevated
## CO2 than at ambient CO2. ## Specifically, noninoculated pots grown under 
## ambient CO2 had 9.5% higher Ncost, while noninoculated pots grown under
## elevated CO2 had 27.6% higher Ncost. Perhaps due to higher Ndemand?

## Individual effect of CO2
emmeans(wpn, pairwise~co2)
## elevated co2 generally has higher Ncost

## Individual effect of inoc
emmeans(wpn, pairwise~inoc)
## Noninoculated pots generally had higher Ncost

## Individual effect of soil N
test(emtrends(wpn, ~1, "n.trt", regrid = "response"))
emmeans(wpn, ~1, "n.trt", at = list(n.trt = 0), type = "response")
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
cld(emmeans(narea, pairwise~inoc*co2))
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
df$chl.mmolm2[25] <- NA
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
df$jmax25.vcmax25[108] <- NA

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
df$gsw[77] <- NA

gsw <- lmer(sqrt(gsw) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
## Only detectable pattern is a one-way negative effect of N fertilization 
## on Ci:Ca:
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
df$pnue[49] <- NA

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

##########################################################################
## Proportion of N in structure
##########################################################################
df$p.structure[c(39, 49, 108, 109, 111)] <- NA

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
## Proportion of N in Rubisco
##########################################################################
df$p.rubisco[c(49)] <- NA

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
## Two-way interaction between inoculation and n.trt
test(emtrends(p.rub, pairwise~inoc, "n.trt"))
## No effect of n.trt in uninoculated pots, negative effect in inoculated
## pots

## Two-way interaction between inoculation and co2
cld(emmeans(p.rub, pairwise~co2*inoc))
## Uninoculated pots have higher p.rub in ambient CO2 treatment
## No diff between CO2 treatments for inoculated pots

## Individual effect of co2
emmeans(p.rub, pairwise~co2)
## Elevated CO2 generally has higher p.rub

## Individual effect of inoc
emmeans(p.rub, pairwise~inoc)
## Inoculated pots generally have higher p.photo

## Individual effect of n.trt
test(emtrends(p.rub, ~1, "n.trt"))
## Increasing soil N fertilization decreases p.rub

##########################################################################
## Proportion of N in bioenergetics
##########################################################################
df$p.bioe[c(49)] <- NA

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
## Two-way interaction between inoculation and n.trt
test(emtrends(p.bioe, pairwise~inoc, "n.trt"))
## No effect of n.trt in uninoculated pots, negative effect in inoculated
## pots

## Individual effect of co2
emmeans(p.bioe, pairwise~co2)
## Elevated CO2 generally has higher p.bioe

## Individual effect of inoc
emmeans(p.bioe, pairwise~inoc)
## Inoculated pots generally have higher p.bioe

## Individual effect of n.trt
test(emtrends(p.bioe, ~1, "n.trt"))
## Increasing soil N fertilization decreases p.bioe

##########################################################################
## Proportion of N in light harvesting
##########################################################################
df$p.lightharv[c(39, 41, 49)] <- NA

p.light <- lmer(sqrt(p.lightharv) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
## Two-way interaction between inoculation and n.trt
test(emtrends(p.light, pairwise~inoc, "n.trt"))
## No effect of n.trt in uninoculated pots, negative effect in inoculated
## pots

## Individual effect of co2
emmeans(p.light, pairwise~co2)
## Elevated CO2 generally has higher p.bioe

## Individual effect of inoc
emmeans(p.light, pairwise~inoc)
## Inoculated pots generally have higher p.bioe

## Individual effect of n.trt
test(emtrends(p.light, ~1, "n.trt"))
## Increasing soil N fertilization decreases p.bioe

##########################################################################
## Table 1: Whole plant traits
##########################################################################
ncost.coefs <- data.frame(summary(ncost)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.ncost = format(Estimate, scientific = TRUE, digits = 3),
         se.ncost = round(Std..Error, digits = 3),
         t.value.ncost = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.ncost, se.ncost, t.value.ncost) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.ncost) %>%
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
         coef.cbg = format(Estimate, scientific = TRUE, digits = 3),
         se.cbg = round(Std..Error, digits = 3),
         t.value.cbg = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.cbg, se.cbg, t.value.cbg) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.cbg) %>%
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
         coef.wpn = format(Estimate, scientific = TRUE, digits = 3),
         se.wpn = round(Std..Error, digits = 3),
         t.value.wpn = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.wpn, se.wpn, t.value.wpn) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.wpn) %>%
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
         coef.tla = format(Estimate, scientific = TRUE, digits = 3),
         se.tla = round(Std..Error, digits = 3),
         t.value.tla = round(t.value, digits = 3)) %>%
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
         coef.tbio = format(Estimate, scientific = TRUE, digits = 3),
         se.tbio = round(Std..Error, digits = 3),
         t.value.tbio = round(t.value, digits = 3)) %>%
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

table1 <- ncost.table %>% full_join(cbg.table) %>% full_join(wpn.table) %>%
  full_join(tla.table) %>% full_join(tbio.table)
write.csv(table1, file = "../working_drafts/tables/NxCO2xI_table1_WP.csv",
         row.names = FALSE)

##########################################################################
## Table 2: Leaf N content
##########################################################################
narea.coefs <- data.frame(summary(narea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.narea = format(Estimate, scientific = TRUE, digits = 3),
         se.narea = round(Std..Error, digits = 3),
         t.value.narea = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.narea, se.narea, t.value.narea) %>%
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
         coef.nmass = format(Estimate, scientific = TRUE, digits = 3),
         se.nmass = round(Std..Error, digits = 3),
         t.value.nmass = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.nmass, se.nmass, t.value.nmass) %>%
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
         coef.marea = format(Estimate, scientific = TRUE, digits = 3),
         se.marea = round(Std..Error, digits = 3),
         t.value.marea = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.marea, se.marea, t.value.marea) %>%
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

table2 <- narea.table %>% full_join(nmass.table) %>% full_join(marea.table)
write.csv(table2, file = "../working_drafts/tables/NxCO2xI_table2_leafN.csv",
          row.names = FALSE)

##########################################################################
## Table 3: Gas exchange
##########################################################################
vcmax.coefs <- data.frame(summary(vcmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.vcmax = format(Estimate, scientific = TRUE, digits = 3),
         se.vcmax = round(Std..Error, digits = 3),
         t.value.vcmax = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax, se.vcmax, t.value.vcmax) %>%
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
         coef.jmax = format(Estimate, scientific = TRUE, digits = 3),
         se.jmax = round(Std..Error, digits = 3),
         t.value.jmax = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.jmax, se.jmax, t.value.jmax) %>%
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
         coef.jvmax = format(Estimate, scientific = TRUE, digits = 3),
         se.jvmax = round(Std..Error, digits = 3),
         t.value.jvmax = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.jvmax, se.jvmax, t.value.jvmax) %>%
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

gsw.coefs <- data.frame(summary(gsw)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.gsw = format(Estimate, scientific = TRUE, digits = 3),
         se.gsw = round(Std..Error, digits = 3),
         t.value.gsw = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.gsw, se.gsw, t.value.gsw) %>%
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

cica.coefs <- data.frame(summary(cica)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.cica = format(Estimate, scientific = TRUE, digits = 3),
         se.cica = round(Std..Error, digits = 3),
         t.value.cica = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.cica, se.cica, t.value.cica) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.cica) %>%
  print(., row.names = FALSE)

cica.table <- data.frame(Anova(cica)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(cica.coefs) %>%
  dplyr::select(treatment, df = Df, coef.cica, 
                chisq.cica = Chisq, pval.cica = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.cica:pval.cica, round, digits = 3),
         chisq.cica = ifelse(chisq.cica < 0.001 & chisq.cica >= 0, 
                            "<0.001", chisq.cica),
         pval.cica = ifelse(pval.cica < 0.001 & pval.cica >= 0, 
                           "<0.001", pval.cica)) %>%
  arrange(treatment)

l.coefs <- data.frame(summary(stomlim)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.l = format(Estimate, scientific = TRUE, digits = 3),
         se.l = round(Std..Error, digits = 3),
         t.value.l = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.l, se.l, t.value.l) %>%
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

table3 <- vcmax.table %>% full_join(jmax.table) %>% full_join(jvmax.table) %>%
  full_join(gsw.table) %>% full_join(cica.table) %>% full_join(l.table)
write.csv(table3, file = "../working_drafts/tables/NxCO2xI_table3_gasEx.csv",
          row.names = FALSE)

##########################################################################
## Table 4: Prop leaf N to photosynthesis, structure, etc.
##########################################################################
p.rub.coefs <- data.frame(summary(p.rub)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.p.rub = format(Estimate, scientific = TRUE, digits = 3),
         se.p.rub = round(Std..Error, digits = 3),
         t.value.p.rub = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.p.rub, se.p.rub, t.value.p.rub) %>%
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
         coef.p.bioe = format(Estimate, scientific = TRUE, digits = 3),
         se.p.bioe = round(Std..Error, digits = 3),
         t.value.p.bioe = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.p.bioe, se.p.bioe, t.value.p.bioe) %>%
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
         coef.p.light = format(Estimate, scientific = TRUE, digits = 3),
         se.p.light = round(Std..Error, digits = 3),
         t.value.p.light= round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.p.light, se.p.light, t.value.p.light) %>%
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
         coef.p.photo = format(Estimate, scientific = TRUE, digits = 3),
         se.p.photo = round(Std..Error, digits = 3),
         t.value.p.photo = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.p.photo, se.p.photo, t.value.p.photo) %>%
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
         coef.p.str = format(Estimate, scientific = TRUE, digits = 3),
         se.p.str = round(Std..Error, digits = 3),
         t.value.p.str = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.p.str, se.p.str, t.value.p.str) %>%
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

table4 <- p.rub.table %>% full_join(p.bioe.table) %>% 
  full_join(p.light.table) %>% full_join(p.photo.table) %>%
  full_join(p.str.table)

write.csv(table4, file = "../working_drafts/tables/NxCO2xI_table4_propN.csv",
          row.names = FALSE)

##########################################################################
## Table 5: PNUE/iWUE
##########################################################################
pnue.coefs <- data.frame(summary(pnue)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.pnue = format(Estimate, scientific = TRUE, digits = 3),
         se.pnue = round(Std..Error, digits = 3),
         t.value.pnue = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.pnue, se.pnue, t.value.pnue) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.pnue) %>%
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
         coef.iwue = format(Estimate, scientific = TRUE, digits = 3),
         se.iwue = round(Std..Error, digits = 3),
         t.value.iwue = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.iwue, se.iwue, t.value.iwue) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.iwue) %>%
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
         coef.narea.gs = format(Estimate, scientific = TRUE, digits = 3),
         se.narea.gs = round(Std..Error, digits = 3),
         t.value.narea.gs = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.narea.gs, se.narea.gs, t.value.narea.gs) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.narea.gs) %>%
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
         coef.vcmax.gs = format(Estimate, scientific = TRUE, digits = 3),
         se.vcmax.gs = round(Std..Error, digits = 3),
         t.value.vcmax.gs = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax.gs, se.vcmax.gs, t.value.vcmax.gs) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.vcmax.gs) %>%
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

table5 <- pnue.table %>% full_join(iwue.table) %>% full_join(narea.gs.table) %>%
  full_join(vcmax.gs.table)
write.csv(table5, file = "../working_drafts/tables/NxCO2xI_table5_pnue_iwue.csv",
          row.names = FALSE)
