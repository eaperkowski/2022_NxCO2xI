library(simr)
library(dplyr)
library(lme4)
library(pwr)
library(car)
library(tidyverse)

df <- data.frame(co2 = rep(rep(c(400, 1000), each = 18), each = 2),
                 soil.n = rep(rep(rep(c(0, 35, 70, 105, 140, 210, 270, 350, 630),
                              times = 2), each = 2), times = 4),
                 inoc = rep(rep(c("yes", "no"), times = 4), times = 18),
                 vcmax = rnorm(144, mean = 76, sd = 5))

pwr.chisq.test(w = 0.45, N = 136, df = 36)
pwr.f2.test(u = 7, v = 136, f2 = 0.25)


test.lm <- lm(vcmax ~ factor(co2) * soil.n * inoc, data = df)
Anova(test.lm)


ggplot(df, aes(x = soil.n, y = vcmax, fill = factor(co2), shape = factor(inoc))) +
  geom_point(size = 2, shape = 21) +
  geom_smooth(method = "lm")



