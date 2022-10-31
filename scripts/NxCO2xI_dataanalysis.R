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
         co2 = as.factor(co2))

# Read in compiled data file

##########################################################################
## Leaf nitrogen contenet (Narea)
##########################################################################
df$nmass.focal[c(113, 114, 117)] <- NA

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
# Stronger stimulation in Narea with increasing soil N in noninoculated pots

# Interaction between CO2 and inoc
emmeans(narea, pairwise~inoc*co2)
# Stronger stimulation in Narea in inoculated pots grown under elevated CO2

# Individual effect of co2
emmeans(narea, pairwise~co2)

# Individual effect of inoc
emmeans(narea, pairwise~inoc)

# Individual effect of n.trt
emtrends(narea, ~1, "n.trt")

##########################################################################
## Vcmax25
##########################################################################
df$p.photo[c(51)] <- NA

p.photo <- lm(p.photo ~ co2 * inoc * n.trt, data = df)

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





