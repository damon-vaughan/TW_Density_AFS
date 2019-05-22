library(tidyverse)

weth <- expand.grid(year = seq(1917,2018), month = seq(1,12)) %>% 
  arrange(year)

Pulliam <- read.csv("RawData/pulliam_weather2.csv", header = TRUE) %>% 
  select(DATE, PRCP, TMIN, TMAX, TAVG) %>% 
  separate(DATE, c("year", "month"), sep = "-", convert = TRUE) %>% 
  rename("PRCPpulliam" = PRCP, "TMINpulliam" = TMIN, "TMAXpulliam" = TMAX, "TAVGpulliam" = TAVG) %>% 
  right_join(weth, by = c("year", "month")) 

FTVjoin <- read.csv("RawData/FortValley_weather.csv", header = TRUE) %>% 
  select(DATE, PRCP, TMIN, TMAX, TAVG) %>% 
  separate(DATE, c("year", "month"), sep = "-", convert = TRUE) %>% 
  rename("PRCPftv" = PRCP, "TMINftv" = TMIN, "TMAXftv" = TMAX, "TAVGftv" = TAVG) %>% 
  right_join(Pulliam, by = c("year", "month"))

# Replace Pulliam na's with values from FTV
Pulliam2 <- FTVjoin %>% 
  replace_na(list(PRCPpulliam = -999, TMINpulliam = -999, TMAXpulliam = -999, TAVGpulliam = -999)) %>% 
  mutate(new.PRCP = ifelse(PRCPpulliam == -999, PRCPftv, PRCPpulliam),
         new.TMIN = ifelse(TMINpulliam == -999, TMINftv, TMINpulliam),
         new.TMAX = ifelse(TMAXpulliam == -999, TMAXftv, TMAXpulliam),
         new.TAVG = ifelse(TAVGpulliam == -999, TAVGftv, TAVGpulliam))
length(which(is.na(Pulliam$PRCPpulliam), arr.ind = T)) - length(which(is.na(Pulliam2$new.PRCP), arr.ind = T))
# Replaced 171 values with FTV. These are mostly at beginning and end, which don't really end up being used. The other big section is years 1943 to 1950. After 1950 the Pulliam data is very consistent.

ggplot(Pulliam2) +
  stat_summary(aes(x = year, y = new.PRCP), fun.y = sum, geom = "line", color = "blue")  

# Adjust to water year
clim <- Pulliam2 %>% 
  select(year, month, PRCP = new.PRCP, TMIN = new.TMIN, TMAX = new.TMAX, TAVG = new.TAVG) %>% 
  mutate(grow.yr = ifelse(month >= 10, year + 1, year)) 

# Get rid of all remaining NAs- Only 51 of them, and mostly at start/end
which(is.na(clim), arr.ind = T)

monthlymean <- clim %>% 
  group_by(month) %>% 
  summarise(PRCP.month = mean(PRCP, na.rm = TRUE),
            TMIN.month = mean(TMIN, na.rm = TRUE),
            TMAX.month = mean(TMAX, na.rm = TRUE),
            TAVG.month = mean(TAVG, na.rm = TRUE))

clim2 <- clim %>%
  left_join(monthlymean, by = "month")

clim2$PRCP[is.na(clim2$PRCP)] <- 
  clim2$PRCP.month[is.na(clim2$PRCP)]
clim2$TMIN[is.na(clim2$TMIN)] <- 
  clim2$TMIN.month[is.na(clim2$TMIN)]
clim2$TMAX[is.na(clim2$TMAX)] <- 
  clim2$TMAX.month[is.na(clim2$TMAX)]
clim2$TAVG[is.na(clim2$TAVG)] <- 
  clim2$TAVG.month[is.na(clim2$TAVG)]

clim3 <- clim2 %>% 
  select(-PRCP.month, -TMIN.month, -TMAX.month, -TAVG.month)

ggplot(climWide) +
  geom_line(aes(x = grow.yr, y = PRCPtot))

# PDSI
PDSI <- read.csv("RawData/PDSI_formatted.csv", header = TRUE) %>% 
  select(YearMonth, PDSI) %>% 
  separate(YearMonth, into = c("year", "month"), sep = 4, convert = TRUE)

# Putting it together
climData <- clim3 %>% 
  left_join(PDSI, by = c("year", "month")) %>% 
  filter(grow.yr >= 1918 & year < 2018) %>% 
  mutate(quarter = ifelse(month >= 10 & month <= 12, 1,
                          ifelse(month >= 1 & month <= 3, 2,
                                 ifelse(month >= 4 & month <= 6, 3,
                                        4))))

write.csv(climData, "AnalysisReadyData/climLong_TWPaper.csv", row.names = FALSE)

climSum <- climData %>% 
  group_by(grow.yr, quarter) %>% 
  summarise(PRCP = sum(PRCP), 
            TMIN = mean(TMIN),
            TMAX = mean(TMAX),
            TAVG = mean(TAVG),
            PDSI = mean(PDSI)) %>%
  filter(grow.yr <= 2017 & grow.yr > 1919) 

climWide <- climSum %>% 
  gather("temporary", "value", PRCP:PDSI) %>% 
  unite(temporary.quarter, c(temporary,quarter), sep = ".") %>% 
  spread(temporary.quarter, value) %>% 
  mutate(PRCPtot = PRCP.1 + PRCP.2 + PRCP.3 + PRCP.4,
         TMINavg = (TMIN.1 + TMIN.2 + TMIN.3 + TMIN.4)/4,
         TMAXavg = (TMAX.1 + TMAX.2 + TMAX.3 + TMAX.4)/4,
         TAVGavg = (TAVG.1 + TAVG.2 + TAVG.3 + TAVG.4)/4,
         PDSIavg = (PDSI.1 + PDSI.2 + PDSI.3 + PDSI.4)/4)

write.csv(climWide, "AnalysisReadyData/climWide_TWPaper.csv", row.names = FALSE)
