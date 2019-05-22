library(tidyverse)

# Weather data used is Global Summaries of the Month from NCEI. Precip is in mm, PRCP includes all precip including snow. Therefore impossible to have zero precip and some snow, except through instrument error. Temp is in Celsius

# Note that there seem to be instrument problems in years 1997-2003 as zero snow is recorded in these datasets. Only worth fixing if I do anything with snow, but I have some additional data downloaded to work with if needed.

# Create dataframe

weth <- expand.grid(year = seq(1917,2018), month = seq(1,12)) %>% 
  arrange(year)

FTV <- read.csv("RawData/FortValley_weather.csv", header = TRUE) %>% 
  select(DATE, PRCP, SNOW, TAVG) %>% 
  separate(DATE, c("year", "month"), sep = "-", convert = TRUE) %>% 
  rename(FTV.prec = PRCP, FTV.snow = SNOW, FTV.temp = TAVG) %>% 
  right_join(weth, by = c("year", "month"))

Pulliam <- read.csv("RawData/pulliam_weather2.csv", header = TRUE) %>% 
  select(DATE, PRCP, SNOW, TAVG) %>% 
  separate(DATE, c("year", "month"), sep = "-", convert = TRUE) %>% 
  rename(pull.prec = PRCP, pull.snow = SNOW, pull.temp = TAVG) %>% 
  right_join(FTV, by = c("year", "month")) 

# Check diffs between the two. The only really concerning one with FTV is 2001
ggplot(Pulliam) +
  stat_summary(aes(x = year, y = pull.prec), fun.y = mean, geom = "line", color = "red") +
  stat_summary(aes(x = year, y = FTV.prec), fun.y = mean, geom = "line", color = "blue")

# Replace FTV na's with values from Pulliam
Pulliam2 <- Pulliam %>% 
  replace_na(list(FTV.prec = -999, FTV.snow = -999)) %>% 
  mutate(new.prec = ifelse(FTV.prec == -999, pull.prec, FTV.prec),
         new.snow = ifelse(FTV.snow == -999, pull.snow, FTV.snow),
         new.temp = ifelse(FTV.temp == -999, pull.temp, FTV.temp))

# Replace snow data from after 1996 with Pulliam's- all zeros in FTV
Pulliam3 <- Pulliam2 %>%
  mutate(new.snow = ifelse(year >= 1997, pull.snow, new.snow))

# Jan 2001 to March 2003: All zero or NA in FTV, suspect instrument problem. Replace 2001 to 2003 with Pulliam
Pulliam4 <- Pulliam3 %>% 
  mutate(new.prec = ifelse(year >= 2001 & year <= 2003, pull.prec, new.prec))

ggplot(Pulliam4) +
  stat_summary(aes(x = year, y = pull.prec), fun.y = sum, geom = "line", color = "red") +
  stat_summary(aes(x = year, y = new.prec), fun.y = sum, geom = "line", color = "blue")  

# Adjust to water year rather than calendar year
clim <- Pulliam4 %>% 
  select(year, month, prec = new.prec, snow = new.snow, temp = new.temp) %>% 
  mutate(grow.yr = ifelse(month >= 10, year + 1, year)) 

# Get rid of all remaining NAs
monthlymean <- clim %>% 
  group_by(month) %>% 
  summarise(prec.month = mean(prec, na.rm = TRUE),
            snow.month = mean(snow, na.rm = TRUE),
            temp.month = mean(temp, na.rm = TRUE))

clim2 <- clim %>%
  left_join(monthlymean, by = "month")

clim2$prec[is.na(clim2$prec)] <- 
  clim2$prec.month[is.na(clim2$prec)]
clim2$snow[is.na(clim2$snow)] <- 
  clim2$snow.month[is.na(clim2$snow)]
clim2$temp[is.na(clim2$temp)] <- 
  clim2$temp.month[is.na(clim2$temp)]

prec.snow.temp <- clim2 %>% 
  select(-prec.month, -snow.month, -temp.month)

ggplot(prec.snow.temp) +
  geom_smooth(aes(x = month, y = temp))

# PDSI --------------------------------------------------------------------

# No NAs but it is missing Dec 2018, so just need to add one row

# Available from: https://www7.ncdc.noaa.gov/CDO/CDODivisionalSelect.jsp. You can find it by searching NCDC divisional data. I'd like to understand more about this dataset. How is PDSI calculated, for one? 

# Updated 12/4/18. Used Arizona Division 2 (northeast) data- this is a huge division, so Flagstaff is in NE.

PDSI <- read.csv("RawData/PDSI_formatted.csv", header = TRUE) %>% 
  select(YearMonth, PDSI) %>% 
  separate(YearMonth, into = c("year", "month"), sep = 4, convert = TRUE)

# Growing Degree Days -----------------------------------------------------

# Base temp is 3.9C. This comes from Nitschke and Innes 2008, who cite the source it came from. This is likely for NW ponderosa, but whatever. GDD is calculated as mean daily temp minus the base temp, then summed. Only days after April 1 are accepted to avoid unseasonably warm winter days.

# This is all to make a blank template for daily climate data

years <- data.frame(seq(1923,2018)) %>% 
  rename("year" = "seq.1923..2018.") %>% 
  mutate(leap = ifelse(year%%4 == 0, "DoMl", "DoM"))

DoM <- data.frame(cbind(
  (c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4, 30), rep(5, 31), rep(6, 30), 
   rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31))),  
  (c(seq(1,31), seq(1,28), seq(1,31), seq(1,30), seq(1,31), seq(1,30), seq(1,31), seq(1,31), 
     seq(1,30), seq(1,31), seq(1,30), seq(1,31))))) 

DoM2 <- purrr::map_df(seq_len(72), ~DoM) %>% 
  rename("month" = "X1", "day" = "X2")

DoMl <- data.frame(cbind(
  (c(rep(1, 31), rep(2, 29), rep(3, 31), rep(4, 30), rep(5, 31), rep(6, 30), 
     rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31))),  
  (c(seq(1,31), seq(1,29), seq(1,31), seq(1,30), seq(1,31), seq(1,30), seq(1,31), seq(1,31), 
     seq(1,30), seq(1,31), seq(1,30), seq(1,31))))) 

DoMl2 <- purrr::map_df(seq_len(24), ~DoMl) %>% 
  rename("month" = "X1", "day" = "X2")

nonleap <- filter(years, leap == "DoM")[,1]
leapy <- filter(years, leap == "DoMl")[,1]

nonleap.df <- data.frame(rep(nonleap, 365)) %>% 
  rename("year" = "rep.nonleap..365.") %>% 
  arrange(year) %>% 
  mutate(month = DoM2$month,
         day = DoM2$day)

leap.df <- data.frame(rep(leapy, 366)) %>% 
  rename("year" = "rep.leapy..366.") %>% 
  arrange(year) %>% 
  mutate(month = DoMl2$month,
         day = DoMl2$day)

Cgrid <- rbind(nonleap.df, leap.df) 

Cgrid <- arrange(Cgrid, year)

# Why Pulliam? FTV would be better, but cant check now cuz of shutdown

PullDaily <- read.csv("RawData/Pulliam_Daily.csv", header = TRUE) %>% 
  mutate(TMAX.c = (TMAX - 32) * (5/9),
         TMIN.c = (TMIN - 32) * (5/9)) %>% 
  select(DATE, TMAX.c, TMIN.c) %>% 
  separate(DATE, c("year", "month", "day"), sep = "-", convert = TRUE) 

PD2 <- PullDaily %>%  
  mutate(DD = (TMAX.c + TMIN.c)/2 - 3.9,
         DD = ifelse(DD < 0, 0, DD),
         DD = ifelse(month < 4, 0, DD)) %>% 
  select(-TMAX.c, -TMIN.c)

PD3 <- left_join(Cgrid, PD2, by = c("year", "month", "day"))

# Replace NAs with daily mean
dailymean <- PD3 %>% 
  group_by(month, day) %>% 
  summarise(DD.day = mean(DD, na.rm = TRUE))

PD4 <- PD3 %>%
  left_join(dailymean, by = c("month", "day"))

PD4$DD[is.na(PD4$DD)] <- 
  PD4$DD.day[is.na(PD4$DD)]

GDDmonth <- PD4 %>% 
  group_by(year, month) %>% 
  summarise(GDDmonth = sum(DD))

GDDyear <- PD4 %>% 
  group_by(year) %>% 
  summarise(GDDyear = sum(DD))

ggplot(GDDyear) +
  geom_point(aes(x = year, y = GDDyear))

# Putting it all together... make final DF's long, summary, and wide --------

# Current df (Jan 3 2019) includes precip, month avg temp, and GDD. Not including PDSI or snow, but these could easily be added.

climData <- prec.snow.temp %>% 
  left_join(GDDmonth, by = c("year", "month")) %>% 
  left_join(PDSI, by = c("year", "month")) %>% 
  select(-snow) %>% 
  filter(grow.yr >= 1918 & year < 2018) %>% 
  mutate(quarter = ifelse(month >= 10 & month <= 12, 1,
                          ifelse(month >= 1 & month <= 3, 2,
                                 ifelse(month >= 4 & month <= 6, 3,
                                        4))))

write.csv(climData, "AnalysisReadyData/climData.csv", row.names = FALSE)
climData <- read.csv("AnalysisReadyData/climData.csv", header = TRUE)

climSum <- climData %>% 
  group_by(grow.yr, quarter) %>% 
  summarise(prec = sum(prec, na.rm = TRUE), 
            temp = mean(temp, na.rm = TRUE),
            GDD = sum(GDDmonth),
            PDSI = mean(PDSI, na.rm = TRUE)) %>%
  filter(grow.yr <= 2017 & grow.yr > 1919) 

climWide <- climSum %>% 
  gather("temporary", "value", prec:PDSI) %>% 
  unite(temporary_quarter, c(temporary,quarter)) %>% 
  spread(temporary_quarter, value) %>% 
  mutate(prectot = prec_1+prec_2+prec_3+prec_4,
         tempavg = (temp_1+temp_2+temp_3+temp_4)/4,
         GDDtot = sum(GDD_1, GDD_2, GDD_3, GDD_4),
         PDSIavg = (PDSI_1 + PDSI_2 + PDSI_3 + PDSI_4)/4)

precJuly <- climData %>% 
  select(grow.yr, month, prec) %>% 
  filter(month == 7) %>% 
  rename("precJuly" = "prec")

write.csv(climWide, "AnalysisReadyData/climWide.csv", row.names = FALSE)
write.csv(precJuly, "AnalysisReadyData/precJuly.csv", row.names = FALSE)

climWide <- read.csv("AnalysisReadyData/climWide.csv", header = TRUE)
climLong <- read.csv("AnalysisReadyData/climData.csv", header = TRUE) %>% 
  filter(grow.yr >= 1919 & grow.yr <= 2017)

# Find mean total values
climSum <- climLong %>%
  group_by(grow.yr) %>%
  summarise(prectot = sum(prec, na.rm = TRUE),
            tempmean = mean(temp, na.rm = TRUE))

# Precip and temp plot
p1 <- ggplot(climSum) +
  geom_line(aes(x = grow.yr, y = prectot)) +
  geom_hline(yintercept = mean(climSum$prectot)) +
  ylab("Total Precipitation (mm)") +
  xlab("Growth Year") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) 

# Summaries
ggplot(climWide) +
  geom_line(aes(x = grow.yr, y = scale(prectot)), color = "red") +
  geom_line(aes(x = grow.yr, y = scale(tempavg)), color = "blue") +
  geom_line(aes(x = grow.yr, y = scale(GDDtot)), color = "green")
  geom_vline(xintercept = 2002) +
  geom_vline(xintercept = 1996) +
  geom_vline(xintercept = 1989) +
  geom_vline(xintercept = 1977) +
  geom_vline(xintercept = 1974) +
  geom_vline(xintercept = 1971) +
  geom_vline(xintercept = 1963) +
  geom_vline(xintercept = 1956) +
  geom_vline(xintercept = 1951)

check <- climData %>% 
  group_by(year) %>% 
  summarise(count = n_distinct(month))
