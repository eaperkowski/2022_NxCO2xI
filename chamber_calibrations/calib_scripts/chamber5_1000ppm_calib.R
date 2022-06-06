###############################################################################
## Calibration script for eCO2 experiment: chamber 5 @ 1000 ppm CO2
###############################################################################
# Chamber program:
#
# Constants:
#   - 50 % relative humidity
#   - 1000 ppm CO2
#
# Program timestep
#   - 10:00   25 deg C;  75 % light intensity
#   - 10:45   21 deg C;  50 % light intensity
#   - 11:30   21 deg C;  25 % light intensity
#   - 12:15   17 deg C;   0 % light intensity
#   - 13:30   25 deg C; 100 % light intensity
#
# Program offsets (before correction in this script)
#   - Day temp point 1:       21.0        Day temp offset 1:        4.0
#   - Day temp point 2:       25.0        Day temp offset 2:        4.0
#   - Day temp point 3:       45.0        Day temp offset 3:        0.0
#
#   - Night temp point 1:     17.0        Night temp offset 1:      4.0
#   - Night temp point 2:     35.0        Night temp offset 2:      0.0 
#   - Night temp point 3:     45.0        Night temp offset 3:      0.0
#
#   - Day humidity offset:     0.0        Night humidity offset:    0.0
#   - Day auxillary offset:  239.0        Night auxillary offset: 200.0

###############################################################################
## Read libraries
###############################################################################
library(readLicorData)
library(tidyverse)
library(dplyr)
library(lubridate)
library(patchwork)
library(ggpubr)

###############################################################################
## Load licor and chamber files
###############################################################################
licor <- licorData("../chamber5_1000ppm_licor") %>%
  mutate(date = ymd(str_match(string = date, 
                              pattern = "[0-9]{4}-[0-9]{2}-[0-9]{2}"))) %>%
  unite(col = "date", date:hhmmss..5, sep = " ") %>%
  mutate(date = strptime(as.POSIXct(date), format = "%Y-%m-%d %H:%M:%S", 
                         tz = "America/Chicago")) %>%
  filter(date > "2022-06-03 16:30:00" & date < "2022-06-04 00:15:00") %>%
  select(date, id, machine, CO2_r, CO2_s, Tair, Tleaf, 
         Txchg, TleafEB, RHcham, Qamb_out)

chamber5 <- read.csv("../chamber5_1000ppm.csv") %>%
  mutate(Day = str_pad(Day, width = 2, pad = "0"),
         Month = str_pad(Month, width = 2, pad = "0"),
         Hour = str_pad(Hour, width = 2, pad = "0"),
         Minute = str_pad(Minute, width = 2, pad = "0")) %>%
  unite(col = "date", Year:Day, sep = "") %>%
  mutate(date = ymd(date)) %>%
  unite(col = "time", Hour:Second, sep = ":") %>%
  unite(col = "date", date:time, sep = " ") %>%
  mutate(date = strptime(as.POSIXct(date), format = "%Y-%m-%d %H:%M:%S", 
                         tz = "America/Chicago")) %>%
  select(date, measured_temp = PV_1, prog_temp = SP_1, 
         measured_rh = PV_2, prog_rh = SP_2,
         measured_co2 = PV_3, prog_co2 = SP_3) %>% 
  filter(date > "2022-06-03 16:30:00" & date < "2022-06-04 00:15:00")

###############################################################################
# Visualize temperature differences between set program, sensor measurements,
# and Licor readings
###############################################################################
# Temperature
temp <- ggplot() +
  geom_line(data = chamber5, aes(x = as.POSIXct(date), 
                                 y = measured_temp, 
                                 color = "measured"),
            size = 0.5) +
   geom_line(data = chamber5, aes(x = as.POSIXct(date), 
                                  y = prog_temp, 
                                  color = "prog"),
             size = 0.5) +
  geom_line(data = licor, aes(x = as.POSIXct(date), 
                              y = as.numeric(Tair), 
                              color = "licor"),
            size = 0.5) +
  scale_color_manual(values = c("red", "blue", "black"),
                     labels = c("Licor",
                                "Chamber sensor", 
                                "Chamber set point")) +
  scale_x_datetime(limits = c(as.POSIXct("2022-06-03 16:30"),
                              as.POSIXct("2022-06-04 00:30")),
                   breaks = c(as.POSIXct("2022-06-03 16:30"),
                              as.POSIXct("2022-06-03 18:30"),
                              as.POSIXct("2022-06-03 20:30"),
                              as.POSIXct("2022-06-03 22:30"),
                              as.POSIXct("2022-06-04 00:30")),
                   date_labels = "%R") +
  labs(x = NULL, y = expression("Air temperature ("~degree~"C)"),
       color = "Measurement type") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        legend.text.align = 0)
temp  

# CO2
co2 <- ggplot() +
  geom_line(data = chamber5, aes(x = as.POSIXct(date), 
                                 y = measured_co2, 
                                 color = "measured"),
            size = 0.5) +
  geom_line(data = chamber5, aes(x = as.POSIXct(date), 
                                 y = prog_co2, 
                                 color = "prog"),
            size = 0.5) +
  geom_line(data = licor, aes(x = as.POSIXct(date), 
                              y = as.numeric(CO2_r), 
                              color = "licor"),
            size = 0.5) +
  scale_color_manual(values = c("red", "blue", "black"),
                     labels = c("Licor",
                                "Chamber sensor", 
                                "Chamber set point")) +
  scale_x_datetime(limits = c(as.POSIXct("2022-06-03 16:30"),
                              as.POSIXct("2022-06-04 00:30")),
                   breaks = c(as.POSIXct("2022-06-03 16:30"),
                              as.POSIXct("2022-06-03 18:30"),
                              as.POSIXct("2022-06-03 20:30"),
                              as.POSIXct("2022-06-03 22:30"),
                              as.POSIXct("2022-06-04 00:30")),
                   date_labels = "%R") +
  scale_y_continuous(limits = c(800, 1200), breaks = seq(800, 1200, 200)) +
  labs(x = "Date", y = expression("CO"[2]~ "(μmol mol"^-1~")"),
       color = "Measurement type") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
co2

# Combine both temp and CO2 plots
temp / co2 + plot_layout(guides = "collect")

###############################################################################
# Day CO2 offsets
###############################################################################
day.licor <- subset(licor, Qamb_out > 1)
day.ch5 <- subset(chamber5,  prog_temp > 17)

day.co2 <- ggplot() +
  geom_density(data = day.ch5, 
               aes(x = as.numeric(measured_co2), fill = "chamber"), alpha = 0.75) +
  geom_density(data = day.licor, 
               aes(x = as.numeric(CO2_r), fill = "licor"), alpha = 0.75) +
  geom_vline(xintercept = 1000, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(700, 1100), breaks = seq(700, 1100, 100)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = expression("CO"[2]~ "(μmol mol"^-1~"CO"[2]~")"),
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
day.co2

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.dayco2.summary <- day.licor %>%
  summarize(co2.mean = mean(as.numeric(CO2_r), na.rm = TRUE),
            co2.ci = 1.96 + (sd(CO2_r)/sqrt(length(CO2_r))),
            co2.uci = co2.mean + co2.ci,
            co2.lci = co2.mean - co2.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.dayco2.summary <- day.ch5 %>%
  summarize(co2.mean = mean(as.numeric(measured_co2), na.rm = TRUE),
            co2.ci = 1.96 + (sd(measured_co2)/sqrt(length(measured_co2))),
            co2.uci = co2.mean + co2.ci,
            co2.lci = co2.mean - co2.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())

# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
dayco2.summary <- li.dayco2.summary %>%
  full_join(ch5.dayco2.summary) %>%
  mutate(co2.offset = (co2.mean[2] - co2.mean[1]),
         co2.offset.actual = 239 - co2.offset) %>%
  data.frame()
dayco2.summary

# Chamber 4 day CO2 offset: 90.7 ppm CO2

###############################################################################
# Night CO2 offsets
###############################################################################
night.licor <- subset(licor, Qamb_out < 1)
night.ch5 <- subset(chamber5, prog_temp == 17)

night.co2 <- ggplot() +
  geom_density(data = night.ch5, 
               aes(x = as.numeric(measured_co2), fill = "chamber"), alpha = 0.75) +
  geom_density(data = night.licor, 
               aes(x = as.numeric(CO2_r), fill = "licor"), alpha = 0.75) +
  geom_vline(xintercept = 1000, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(700, 1100), breaks = seq(700, 1100, 100)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = expression("CO"[2]~ "(μmol mol"^-1~"CO"[2]~")"),
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
night.co2

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.nightco2.summary <- night.licor %>%
  summarize(co2.mean = mean(as.numeric(CO2_r), na.rm = TRUE),
            co2.ci = 1.96 + (sd(CO2_r)/sqrt(length(CO2_r))),
            co2.uci = co2.mean + co2.ci,
            co2.lci = co2.mean - co2.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.nightco2.summary <- night.ch5 %>%
  summarize(co2.mean = mean(as.numeric(measured_co2), na.rm = TRUE),
            co2.ci = 1.96 + (sd(measured_co2)/sqrt(length(measured_co2))),
            co2.uci = co2.mean + co2.ci,
            co2.lci = co2.mean - co2.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())

# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
nightco2.summary <- li.nightco2.summary %>%
  full_join(ch5.nightco2.summary) %>%
  mutate(co2.offset = (co2.mean[2] - co2.mean[1]),
         co2.offset.actual = 200 - co2.offset) %>%
  data.frame()
nightco2.summary

# Chamber 4 night CO2 offset: 58.9 ppm CO2

###############################################################################
# Day RH offsets
###############################################################################
day.rh <- ggplot() +
  geom_density(data = day.ch5, 
               aes(x = as.numeric(measured_rh), fill = "chamber"), alpha = 0.75) +
  geom_density(data = day.licor, 
               aes(x = as.numeric(RHcham), fill = "licor"), alpha = 0.75) +
  geom_vline(xintercept = 50, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(30, 60), breaks = seq(30, 60, 10)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = "Relative humidity (%)",
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
day.rh

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.dayrh.summary <- day.licor %>%
  summarize(rh.mean = mean(as.numeric(RHcham), na.rm = TRUE),
            rh.ci = 1.96 + (sd(RHcham)/sqrt(length(RHcham))),
            rh.uci = rh.mean + rh.ci,
            rh.lci = rh.mean - rh.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())
li.dayrh.summary

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.dayrh.summary <- day.ch5 %>%
  summarize(rh.mean = mean(as.numeric(measured_rh), na.rm = TRUE),
            rh.ci = 1.96 + (sd(measured_rh)/sqrt(length(measured_rh))),
            rh.uci = rh.mean + rh.ci,
            rh.lci = rh.mean - rh.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())

# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
dayrh.summary <- li.dayrh.summary %>%
  full_join(ch5.dayrh.summary) %>%
  mutate(rh.offset = (rh.mean[2] - rh.mean[1]),
         rh.offset.actual = rh.offset) %>%
  data.frame()
dayrh.summary

# Chamber 4 day RH offset: -10.7 %

###############################################################################
# Night RH offsets
###############################################################################
night.rh <- ggplot() +
  geom_density(data = night.ch5, 
               aes(x = as.numeric(measured_rh), fill = "chamber"), alpha = 0.75) +
  geom_density(data = night.licor, 
               aes(x = as.numeric(RHcham), fill = "licor"), alpha = 0.75) +
  geom_vline(xintercept = 50, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(30, 60), breaks = seq(30, 60, 10)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = "Relative humidity (%)",
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
night.rh

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.nightrh.summary <- night.licor %>%
  summarize(rh.mean = mean(as.numeric(RHcham), na.rm = TRUE),
            rh.ci = 1.96 + (sd(RHcham)/sqrt(length(RHcham))),
            rh.uci = rh.mean + rh.ci,
            rh.lci = rh.mean - rh.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())
li.nightrh.summary

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.nightrh.summary <- night.ch5 %>%
  summarize(rh.mean = mean(as.numeric(measured_rh), na.rm = TRUE),
            rh.ci = 1.96 + (sd(measured_rh)/sqrt(length(measured_rh))),
            rh.uci = rh.mean + rh.ci,
            rh.lci = rh.mean - rh.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())

# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
nightrh.summary <- li.nightrh.summary %>%
  full_join(ch5.nightrh.summary) %>%
  mutate(rh.offset = (rh.mean[2] - rh.mean[1]),
         rh.offset.actual = rh.offset) %>%
  data.frame()
nightrh.summary

# Chamber 4 night RH offset: -15.2 %

###############################################################################
# 25 deg C day offsets
###############################################################################
# Subset Licor and chamber measurements that were set at 25degC
li.25C <- subset(licor, date > "2022-06-03 22:30:00")
ch5.25C <- subset(chamber5, date > "2022-06-03 22:30:00")

# Visualize density plots of chamber sensor temperature and licor Tair
ch5.25C$measured_temp

dens.25C <- ggplot() +
  geom_density(data = ch5.25C, 
               aes(x = measured_temp, 
                   fill = "chamber"), alpha = 0.75) +
  geom_density(data = li.25C, 
               aes(x = as.numeric(Tair), fill = "licor"), 
               alpha = 0.75) +
  geom_vline(xintercept = 25, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(24.5, 26), breaks = seq(24.5,26,0.5)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = expression("Air temperature ("~degree~"C)"),
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
dens.25C

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.25C.summary <- li.25C %>%
  summarize(temp.mean = mean(as.numeric(Tair)),
            temp.ci = 1.96 + (sd(Tair)/sqrt(length(Tair))),
            temp.uci = temp.mean + temp.ci,
            temp.lci = temp.mean - temp.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.25C.summary <- ch5.25C %>%
  summarize(temp.mean = mean(as.numeric(measured_temp)),
            temp.ci = 1.96 + (sd(measured_temp)/sqrt(length(measured_temp))),
            temp.uci = temp.mean + temp.ci,
            temp.lci = temp.mean - temp.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())


# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
temp.25C.summary <- li.25C.summary %>%
  full_join(ch5.25C.summary) %>%
  mutate(temp.offset = temp.mean[2] - temp.mean[1],
         temp.offset.actual = 4 - temp.offset) %>%
  data.frame()
temp.25C.summary

###############################################################################
# 21 deg C day offsets
###############################################################################
# Subset Licor and chamber measurements that were set at 25degC
li.21C <- subset(licor, date > "2022-06-03 18:30:00" & date < "2022-06-03 19:15:00")
ch5.21C <- subset(chamber5, date > "2022-06-03 18:30:00" & date < "2022-06-03 19:15:00")

# Visualize density plots of chamber sensor temperature and licor Tair
dens.21C <- ggplot() +
  geom_density(data = ch5.21C, 
               aes(x = as.numeric(measured_temp), 
                   fill = "chamber"), alpha = 0.75) +
  geom_density(data = li.21C, 
               aes(x = as.numeric(Tair), fill = "licor"), 
               alpha = 0.75) +
  geom_vline(xintercept = 21, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(20, 23), breaks = seq(20, 23, 1)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = expression("Air temperature ("~degree~"C)"),
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
dens.21C

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.21C.summary <- li.21C %>%
  summarize(temp.mean = mean(as.numeric(Tair)),
            temp.ci = 1.96 + (sd(Tair)/sqrt(length(Tair))),
            temp.uci = temp.mean + temp.ci,
            temp.lci = temp.mean - temp.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.21C.summary <- ch5.21C %>%
  summarize(temp.mean = mean(as.numeric(measured_temp)),
            temp.ci = 1.96 + (sd(measured_temp)/sqrt(length(measured_temp))),
            temp.uci = temp.mean + temp.ci,
            temp.lci = temp.mean - temp.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())

# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
temp.21C.summary <- li.21C.summary %>%
  full_join(ch5.21C.summary) %>%
  mutate(temp.offset = temp.mean[2] - temp.mean[1],
         temp.offset.actual = 4 - temp.offset) %>%
  data.frame()
temp.21C.summary

###############################################################################
# 17 deg C night offsets
###############################################################################
# Subset Licor and chamber measurements that were set at 25degC
li.17C <- subset(licor, date > "2022-06-03 20:30:00" & date < "2022-06-03 21:30:00")
ch5.17C <- subset(chamber5, date > "2022-06-03 20:30:00" & date < "2022-06-03 21:30:00")

# Visualize density plots of chamber sensor temperature and licor Tair
dens.17C <- ggplot() +
  geom_density(data = ch5.17C, 
               aes(x = as.numeric(measured_temp), 
                   fill = "chamber"), alpha = 0.75) +
  geom_density(data = li.17C, 
               aes(x = as.numeric(Tair), fill = "licor"), 
               alpha = 0.75) +
  geom_vline(xintercept = 17, size = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(15, 25), breaks = seq(15, 25, 2)) +
  scale_fill_brewer(palette = "Spectral", labels = c("Chamber w/o offset",
                                                     "Licor")) +
  labs(x = expression("Air temperature ("~degree~"C)"),
       y = "Density", fill = "Measurement type") +
  theme_bw(base_size = 18)
dens.17C

# Calculate mean, ci, uci, and lci of licor measurements at 25degC
li.17C.summary <- li.17C %>%
  summarize(temp.mean = mean(as.numeric(Tair)),
            temp.ci = 1.96 + (sd(Tair)/sqrt(length(Tair))),
            temp.uci = temp.mean + temp.ci,
            temp.lci = temp.mean - temp.ci) %>%
  mutate(meas.type = "licor") %>%
  select(meas.type, everything())

# Calculate mean, ci, uci, and lci of chamber sensor measurements at 25degC
ch5.17C.summary <- ch5.17C %>%
  summarize(temp.mean = mean(as.numeric(measured_temp)),
            temp.ci = 1.96 + (sd(measured_temp)/sqrt(length(measured_temp))),
            temp.uci = temp.mean + temp.ci,
            temp.lci = temp.mean - temp.ci) %>%
  mutate(meas.type = "chamber.sensor") %>%
  select(meas.type, everything())

# Merge licor and sensor dataframes. Add offsets and correct
# from original offsets before program set (see top of code)
temp.17C.summary <- li.17C.summary %>%
  full_join(ch5.17C.summary) %>%
  mutate(temp.offset = temp.mean[2] - temp.mean[1],
         temp.offset.actual = 5 - temp.offset) %>%
  data.frame()
temp.17C.summary
