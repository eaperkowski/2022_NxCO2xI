##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(multcomp)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Read in compiled data file
df <- read.csv("../data_sheets/NxCO2xI_compiled_datasheet.csv") %>%
  mutate(n.trt = as.numeric(n.trt),
         inoc = factor(inoc, levels = c("no.inoc", "inoc")),
         co2 = factor(co2, levels = c("amb", "elv")),
         co2.inoc = str_c(co2, "_", inoc),
         nod.root.ratio = nodule.biomass / root.biomass,
         pnue.growth = anet.growth / (narea / 14),
         lar = tla / total.biomass,
         lmf = leaf.biomass / total.biomass,
         smf = stem.biomass / total.biomass,
         rmf = root.biomass / total.biomass) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05))

##########################################################################
## Narea - Anet,420
##########################################################################
anet.420 <- lmer(anet ~ co2 * inoc * narea + (1|rack:co2), data = df)

# Check model assumptions
plot(anet.420)
qqnorm(residuals(anet.420))
qqline(residuals(anet.420))
densityPlot(residuals(anet.420))
shapiro.test(residuals(anet.420))
outlierTest(anet.420)

# Model results
summary(anet.420)
Anova(anet.420)
r.squaredGLMM(anet.420)

# Pairwise comparisons
test(emtrends(anet.420, pairwise~co2, "narea"))
test(emtrends(anet.420, pairwise~inoc, "narea"))
cld(emtrends(anet.420, pairwise~co2*inoc, "narea"))

##########################################################################
## Narea - Anet,420
##########################################################################
anet.gc <- lmer(anet.growth ~ co2 * inoc * narea + (1|rack:co2), data = df)

# Check model assumptions
plot(anet.gc)
qqnorm(residuals(anet.gc))
qqline(residuals(anet.gc))
densityPlot(residuals(anet.gc))
shapiro.test(residuals(anet.gc))
outlierTest(anet.gc)

# Model results
summary(anet.gc)
Anova(anet.gc)
r.squaredGLMM(anet.gc)

# Pairwise comparisons
test(emtrends(anet.gc, pairwise~inoc, "narea"))

##########################################################################
## Vcmax25
##########################################################################
vcmax <- lmer(vcmax25 ~ co2 * inoc * narea + (1|rack:co2), data = df)

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
test(emtrends(vcmax, pairwise~inoc, "narea"))
test(emtrends(vcmax, pairwise~co2*inoc, "narea"))




