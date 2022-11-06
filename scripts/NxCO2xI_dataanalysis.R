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
         co2 = ifelse(co2 == "a", "amb", "elv"),
         inoc = ifelse(inoc == "n", "no.inoc", "inoc"))

##########################################################################
## Leaf nitrogen content (Narea)
##########################################################################
narea <- lm(narea ~ co2 * inoc * n.trt, data = df)

# Check model assumptions
qqnorm(residuals(narea))
qqline(residuals(narea))
densityPlot(residuals(narea))
shapiro.test(residuals(narea))
outlierTest(narea)

# Model results
summary(narea)
Anova(narea)

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
marea <- lm(log(marea) ~ co2 * inoc * n.trt, data = df)

# Check model assumptions
qqnorm(residuals(marea))
qqline(residuals(marea))
densityPlot(residuals(marea))
shapiro.test(residuals(marea))
outlierTest(marea)

# Model results
summary(marea)
Anova(marea)

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
nmass <- lm(nmass.focal ~ co2 * inoc * n.trt, data = df)

# Check model assumptions
qqnorm(residuals(nmass))
qqline(residuals(nmass))
densityPlot(residuals(nmass))
shapiro.test(residuals(nmass))
outlierTest(nmass)

# Model results
summary(nmass)
Anova(nmass)

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
chl.area <- lm(chl.mmolm2 ~ co2 * inoc * n.trt, data = df)

# Check model assumptions
qqnorm(residuals(chl.area))
qqline(residuals(chl.area))
densityPlot(residuals(chl.area))
shapiro.test(residuals(chl.area))
outlierTest(chl.area)

# Model results
summary(chl.area)
Anova(chl.area)

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
df$p.photo[c(51)] <- NA

p.photo <- lm(vcmax25 ~ co2 * inoc * n.trt, data = df)

# Check model assumptions
qqnorm(residuals(p.photo))
qqline(residuals(p.photo))
densityPlot(residuals(p.photo))
shapiro.test(residuals(p.photo))
outlierTest(p.photo)

# Model results
summary(p.photo)
Anova(p.photo)

# Pairwise comparisons

# Two-way interaction between inoculation and n.trt
test(emtrends(p.photo, ~inoc, "n.trt"))
## Increasing N fertilization decreased p.photo in inoculated pots; did
## not influence non-inoculated pots

# Individual effect of CO2
emmeans(p.photo, pairwise~co2)

# Individual effect of inoculation
emmeans(p.photo, pairwise~inoc)

# Individual effect of n.trt
test(emtrends(p.photo, ~1, "n.trt"))





