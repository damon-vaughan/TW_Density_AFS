library(tidyverse)
library(dplR)
library(TTR)
library(treeclim)


# Updated GPNA chronology from Pete ---------------------------------------

GPNA <- read.csv("RawData/GPNA_Updated.csv", header = FALSE, stringsAsFactors = FALSE) %>% 
  rename(SampleID = V1, decade = V2) %>% 
  gather(key = "PlusYear", value = "RW", V3:V12) %>% 
  na.omit()
ConvertYear <- data.frame(PlusYear = c("V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"),
                          Add = seq(0,9), stringsAsFactors = FALSE)
GPNA2 <- left_join(GPNA, ConvertYear, by = "PlusYear") %>% 
  mutate(year = decade + Add) %>% 
  arrange(SampleID, year) %>% 
  select(SampleID, year, RW) %>% 
  mutate(RW = RW/1000) 

# If below <5 it is code for NA, just remove. One got entered as negative, so if <0, flip the sign.
GPNA3 <- GPNA2 %>% 
  filter(RW >= -5) %>% 
  mutate(RW = ifelse(RW <= 0, RW * -1, RW))

# There were two years that were duplicated in sample 262-407 (1940 and 1941). I just took the min
pearsonLong <- GPNA3 %>% 
  group_by(SampleID, year) %>% 
  summarise(RW = min(RW))

pearsonFull <- pearsonLong %>% 
  spread(key = SampleID, value = RW)

# Greybill chrono- raw data as rwl file is available from https://www.ncdc.noaa.gov/paleo-search/?dataTypeId=18
pearsonFull <- read.rwl("RawData/az521_Pearson.rwl")
climWide <- read.csv("AnalysisReadyData/climWide.csv", header = TRUE)

pearsonLong <- pearsonFull %>% 
  mutate(year = as.numeric(row.names(pearsonFull))) %>% 
  select(year, everything()) %>% 
  gather("SampleID", "RW", 2:68) %>% 
  na.omit()

# Calculate BAI- assumes ring series go all the way to pith, which may not be true in all cases
pearsonLong2 <- pearsonLong %>% 
  group_by(SampleID) %>% 
  mutate(mmRingStart = cumsum(RW) - RW + 1,
         BAI = pi * ((mmRingStart + RW)^2 - mmRingStart^2)/100)

ggplot(pearsonLong2) +
  geom_smooth(aes(x = year, y = BAI))

# Most trees minYear from 1900-1950. Want to subset these to look more closely at second growth
pearsonSum <- pearsonLong %>% 
  group_by(SampleID) %>% 
  summarise(length = max(year - min(year)),
            minYR = min(year),
            min = min(RW), max = max(RW))

# The long term trend is weird- big spike around 1900. This could reflect the abundance of trees that start in this year, but many individual old growth trees show this as well.
ggplot(pearsonLong) +
  geom_smooth(aes(x = year, y = RW))

ggplot(filter(pearsonLong, SampleID == "GPN011")) +
  geom_smooth(aes(x = year, y = RW))

# Subset to second growth
pearson1900 <- pearsonFull %>% 
  select(c(filter(pearsonSum, minYR >= 1900)$SampleID)) 

pearson1900.long <- pearson1900 %>% 
  mutate(year = as.numeric(row.names(pearson1900))) %>%
  select(year, everything()) %>% 
  gather("SampleID", "RW", 2:25) %>% 
  na.omit() %>% 
  group_by(SampleID) %>% 
  mutate(cambialAge = year - min(year) + 1)

pearson1900.long <- pearson1900 %>% 
  select(year, everything()) %>% 
  gather("SampleID", "RW", 2:516) %>% 
  na.omit() 
%>% 
  group_by(SampleID) %>% 
  mutate(cambialAge = year - min(year) + 1)

# lots of variability within each tree. This is the climate signal
ggplot(pearson1900.long) +
  geom_line(aes(x = year, y = RW, color = SampleID)) +
  theme(legend.position="none")

# overall trend looks like a negative exponential- need to remove this
ggplot(pearson1900.long) +
  geom_smooth(aes(x = year, y = RW)) +
  theme(legend.position="none")  

# Check heteroscedasticity of raw
dflist <- split(pearson1900.long, f = pearson1900.long$SampleID)
movingSD <- lapply(dflist, function(x){
  runSD(x[,3], n = 10)}) %>%
  unlist() %>% 
  data.frame()

movingAvg <- lapply(dflist, function(x){
  runMean(x[,3], n = 10)}) %>% 
  unlist() %>% 
  data.frame()

pearson1900.long2 <- pearson1900.long %>% 
  ungroup() %>% 
  mutate(movingSD = movingSD[,1],
         movingAvg = movingAvg[,1]) %>% 
  na.omit()
cor(pearson1900.long2[,5:6])

ggplot(pearson1900.long2) +
  geom_point(aes(x = movingSD, y = movingAvg))

# Calculate BAI- assumes ring series go all the way to pith, which may not be true in all cases
pearson1900.long2 <- pearson1900.long %>% 
  group_by(SampleID) %>% 
  mutate(mmRingStart = cumsum(RW) - RW + 1,
         BAI = pi * ((mmRingStart + RW)^2 - mmRingStart^2)/100) %>% 
  filter(cambialAge >= 20)

pearsonChron <- pearson1900.long2 %>% 
  group_by(year) %>% 
  summarise(BAI = mean(BAI))

ggplot(pearson1900.long2) +
  geom_smooth(aes(x = year, y = RW))
ggplot(pearsonChron) +
  geom_point(aes(x = year, y = BAI))
# Detrending --------------------------------------------------------------

# Detrend with negative exponential
pearson1900.detrend <- detrend(pearson1900, method = "ModNegExp")

# Detrend with 67 yr cubic spline (no power transformation), as in Helama paper
pearson1900.detrend <- detrend(pearson1900, method = "Spline")
pearsonFull.detrend <- detrend(pearsonFull, method = "Spline")

pearson1900.detrend.long <- pearson1900.detrend %>% 
  mutate(year = as.numeric(row.names(pearson1900.detrend))) %>% 
  select(year, everything()) %>% 
  gather("SampleID", "RWI", 2:25) %>% 
  na.omit()%>% 
  group_by(SampleID) %>% 
  mutate(cambialAge = year - min(year) + 1)

ggplot(pearson1900.detrend.long) +
  # geom_smooth(aes(x = year, y = RWI), method = "gam", 
              # formula = y ~ s(x, k = 6)) +
  geom_point(aes(x = year, y = RWI)) +
  theme(legend.position="none")  

# Check heteroscedasticity of detrended. The 67n spline also removes the correlation between mean and SD, thus stabilizing the variance. Maybe no need for a power transportation?
dflist <- split(pearson1900.detrend.long, f = pearson1900.detrend.long$SampleID)
movingSD <- lapply(dflist, function(x){
  runSD(x[,3], n = 10)}) %>%
  unlist() %>% 
  data.frame()

movingAvg <- lapply(dflist, function(x){
  runMean(x[,3], n = 10)}) %>% 
  unlist() %>% 
  data.frame()

pearson1900.long2 <- pearson1900.detrend.long %>% 
  ungroup() %>% 
  mutate(movingSD = movingSD[,1],
         movingAvg = movingAvg[,1]) %>% 
  na.omit()
cor(pearson1900.long2[,5:6])

ggplot(pearson1900.long2) +
  geom_point(aes(x = movingSD, y = movingAvg))

# Chronologies ------------------------------------------------------------

pearson1900.chron <- chron(pearson1900.detrend) %>% 
  rename("RWI" = "xxxstd")
pearsonFull.chron <- chron(pearsonFull.detrend)
crn.plot(pearson1900.chron)

# This looks good- centered around 1 with lots of high freq variation. But there is a weird pickup at the end?
ggplot(pearson1900.detrend.long) +
  geom_line(aes(x = year, y = RWI, color = SampleID)) +
  theme(legend.position="none")

chron2 <- pearson1900.chron %>% 
  mutate(year = as.numeric(row.names(pearson1900.chron))) 

chron3 <- left_join(chron2, climWide,
                    by = c("year" = "grow.yr")) %>%
  na.omit() 

cor(chron3)
# Pretty good correlation (.45) between this index and total precip. Might work better quarterly. 

mod <- lm(RWI ~ prectot + GDDtot, data = chron3)
mod <- lm(RWI ~ prectot + temp_4, data = chron3)
summary(mod)

# Very interesting... it appears that droughts reduce growth but high precip years don't necessarily have an equivalent increase. So the linear correlation might not be the best way to represent it.
ggplot(chron3) +
  geom_point(aes(x = prectot, y = RWI))

library(mgcv)
mod <- gam(RWI ~ s(prectot) + s(GDDtot), data = chron3)
summary(mod)

mod <- lm(log(RWI) ~ prcptot, data = chron.correlation)
summary(mod)
summary(climWide$prcptot)
summary(climWide$GDDtot)

# Filter to only include prcp below the median
mod <- lm(RWI ~ prcptot, data = filter(chron.correlation, prcptot <= 550))
summary(mod)
# conclusion: Below the median value, a 100 mm reduction in precip causes a .2 reduction in RWI.

ggplot(filter(chron.correlation, prcptot <= 550)) +
  geom_smooth(aes(x = prcptot, y = RWI), method = "lm", se = FALSE) +
  geom_point(aes(x = prcptot, y = RWI))

# Current working model to remove climate from TW trees. Room for improvment
climMod <- lm(RWI ~ prectot + temp_4, data = chron3)
summary(climMod)

# Adjust mean to zero, then fit a no intercept model
chron4 <- chron3 - 1
mean(chron3$RWI) + mean(chron4$RWI)
mean(chron4$RWI)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.0925740  0.4975142   4.206 8.65e-05 ***
# prectot      0.0007214  0.0001670   4.321 5.83e-05 ***
# temp_4      -0.0990936  0.0316145  -3.134  0.00265 ** 
# Multiple R-squared:  0.3307,	Adjusted R-squared:  0.3087 

# Treeclim ----------------------------------------------------------------

# Pre-packaged spain dataset. Using "response", coeffs depend on what else is in the model. Using "correlation", they don't. Overall, coefficients and significance are low.
spai020 <- spai020
spain_prec <- spain_prec
spain_temp <- spain_temp

spaiboth <- dcc(chrono = spai020,
                climate = list(spain_temp, spain_prec),
                var_names = c("temp", "prec"),
                method = "correlation")
spaitemp <- dcc(chrono = spai020,
                climate = list(spain_temp),
                var_names = c("temp"),
                method = "correlation")
plot(spaiboth)
plot(spaitemp)

precip <- read.csv("AnalysisReadyData/climData.csv", header = TRUE) %>% 
  na.omit() %>% 
  select(year, month, prec) %>% 
  spread(month, prec) %>% 
  rename("JAN" = "1", "FEB" = "2", "MAR" = "3", "APR" = "4", "MAY" = "5",
         "JUN" = "6", "JUL" = "7", "AUG" = "8", "SEP" = "9", "OCT" = "10",
         "NOV" = "11", "DEC" = "12") %>% 
  filter(year > 1917 & year < 2017)

temp <- read.csv("AnalysisReadyData/climData.csv", header = TRUE) %>% 
  select(year, month, temp) %>% 
  spread(month, temp) %>% 
  rename("JAN" = "1", "FEB" = "2", "MAR" = "3", "APR" = "4", "MAY" = "5",
         "JUN" = "6", "JUL" = "7", "AUG" = "8", "SEP" = "9", "OCT" = "10",
         "NOV" = "11", "DEC" = "12") %>% 
  filter(year > 1917 & year < 2017)

GDD <- read.csv("AnalysisReadyData/climData.csv", header = TRUE) %>% 
  select(year, month, GDDmonth) %>% 
  spread(month, GDDmonth) %>% 
  rename("JAN" = "1", "FEB" = "2", "MAR" = "3", "APR" = "4", "MAY" = "5",
         "JUN" = "6", "JUL" = "7", "AUG" = "8", "SEP" = "9", "OCT" = "10",
         "NOV" = "11", "DEC" = "12") %>% 
  filter(year > 1917 & year < 2017) %>% 
  na.omit()

pearson.year <- dcc(chrono = pearson1900.chron,
                    climate = list(precip, temp),
                    var_names = c("prec", "temp"),
                    selection = .sum("prec", -10:9) + 
                      .mean("temp", -10:12) + 
                      .mean("temp", 1:3) + 
                      .mean("temp", 4:6) + 
                      .mean("temp", 7:9),
                    boot = "std")
plot(pearson.year)

pearson.year <- dcc(chrono = pearson1900.chron,
                    climate = list(precip, GDD),
                    var_names = c("prec", "temp"),
                    selection = .sum("temp", -10:9) +
                      .range("temp", 5:9),
                    boot = "std")
plot(pearson.year)
  
pearson.month <- dcc(chrono = pearsonFull.chron,
                     climate = list(precip),
                     var_names = "RWI",
                     boot = "std")
pearson.season <- dcc(chrono = pearsonFull.chron,
                    climate = list(precip),
                    var_names = "xxxstd",
                    selection = .sum("xxxstd", 1:3) + .sum("xxxstd", 4:6) +
                      .sum("xxxstd", 7:9) + .sum("xxxstd", -10:-12),
                    boot = "std")
pearson.select <- dcc(chrono = pearson1900.chron,
                     climate = list(precip),
                     var_names = "RWI",
                     boot = "std",
                     selection = .range("RWI", -10:-12) +
                       .range("RWI", 2:3) +.range("RWI", 5:7))
plot(pearson.month)
plot(pearson.season)
plot(pearson.select)

# Pearson with all climate variables

pearson.select <- dcc(chrono = pearson1900.chron,
                      climate = list(precip, GDD),
                      var_names = c("prec", "GDD"),
                      selection = .range("prec", 5:7) + .range("GDD", -10:9),
                      boot = "std")
plot(pearson.select)
