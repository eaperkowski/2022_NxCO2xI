## Load libraries
library(dplyr)
library(lubridate)

## Read in chamber climate data
chamber.data <- list.files(path = "../chamber_climate",
           pattern = "\\.csv$", recursive = TRUE,
           full.names = TRUE)
chamber.data <- setNames(chamber.data, chamber.data)

## Merge chamber climate data to central df, create month column,
## convert to time unit
chamber.clim <- lapply(chamber.data, read.csv) %>%
  reshape::merge_all() %>%
  mutate(month = str_pad(month, 2, side = "left", pad = "0"),
         hour = str_pad(hour, 2, side = "left", pad = "0"),
         minute = str_pad(minute, 2, side = "left", pad = "0")) %>%
  unite(col = "date", year:day, sep = "-") %>%
  unite(col = "time", hour:second, sep = ":") %>%
  unite(col = "date.time", date:time, sep = " ") %>%
  mutate(date = ymd_hms(date.time))


############ 
##Read in HOBO sensor temp data 
############
hobo.data <- list.files(path = "../hobo_data",
                           pattern = "\\.csv$", recursive = TRUE,
                           full.names = TRUE)
hobo.data <- setNames(hobo.data, hobo.data)

hobo.clim <- lapply(hobo.data, read.csv) %>%
  reshape::merge_all() %>%
  mutate(date.time = mdy_hm(date.time)) %>%
  arrange(chamber, date.time) %>%
  separate(date.time, c("date", "time"), remove = FALSE,
           sep = " ") %>%
  mutate(hour = hour(date.time),
         tod = ifelse(hour >= 1 & hour < 9,
                      "night", "day"),
         temp.set = ifelse(tod == "night", 17,
                           ifelse(hour > 11 & hour < 22,
                                  25,
                                  21))) %>%
  dplyr::select(chamber, date.time, tod, air.temp, temp.set)


eco2.hobo <- hobo.clim %>% 
  filter(date.time > "2022-06-18" & date.time < "2022-08-05")  %>%
  mutate(co2.treat = "elevated")

aco2.hobo <- hobo.clim %>% 
  filter(date.time > "2022-08-14" & date.time < "2022-10-01")  %>%
  mutate(co2.treat = "ambient")

ggplot() +
  geom_line(data = eco2, 
            aes(x = date, y = temp.meas),
            color = "red") +
  geom_line(data = eco2.hobo, 
            aes(x = date.time, y = air.temp),
            color = "blue") +
  facet_wrap(~chamber)

#############
## Filter chamber.clim to dates of experiment, tack on 
## whether this was for elevated or ambient CO2 iteration
#############
eco2 <- chamber.clim %>%
  filter(date > "2022-06-18" & date < "2022-08-05") %>%
  mutate(co2.treat = "elevated")
aco2 <- chamber.clim %>%
  filter(date > "2022-08-14" & date < "2022-10-01") %>%
  mutate(co2.treat = "ambient")

trt.night.means <- eco2 %>% full_join(aco2) %>%
  filter(temp.set == 17) %>%
  group_by(co2.treat, chamber) %>%
  summarize(temp = mean(temp.meas),
            rh = mean(humidity.meas)) %>%
  ungroup(chamber) %>%
  summarize(temp.mean = mean(temp),
            temp.sd = sd(temp),
            rh.mean = mean(rh),
            rh.sd = sd(rh))

trt.day.means <- eco2 %>% full_join(aco2) %>%
  filter(temp.set != 17) %>%
  group_by(co2.treat, chamber) %>%
  summarize(temp = mean(temp.meas),
            rh = mean(humidity.meas)) %>%
  ungroup(chamber) %>%
  summarize(temp.mean = mean(temp),
            temp.sd = sd(temp),
            rh.mean = mean(rh),
            rh.sd = sd(rh))

trt.all.means <- eco2 %>% full_join(aco2) %>%
  group_by(co2.treat, chamber) %>%
  summarize(temp = mean(temp.meas),
            rh = mean(humidity.meas)) %>%
  ungroup(chamber) %>%
  summarize(temp.mean = mean(temp),
            temp.sd = sd(temp),
            rh.mean = mean(rh),
            rh.sd = sd(rh))


co2.variability <- eco2 %>% full_join(aco2) %>%
  group_by(co2.treat, chamber) %>%
  summarize(co2.ppm = mean(co2.meas)) %>%
  ungroup(chamber) %>%
  summarize(co2.mean = mean(co2.ppm),
            co2.se = sd(co2.ppm)/sqrt(6))
  
  
  
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
  
  
  
  
  