library(needs)
needs(tidyverse, car, nlme, mgcv, gridExtra, dplR, treeclim, emmeans, dichromat)

source("furnival.R")
source("bias function.R")

# Read in data and summarise --------------------------------------------------------

df1 <- read.csv("AnalysisReadyData/sumProfile_AFS.csv", 
                header = TRUE, stringsAsFactors = FALSE) 

Add58to62 <- df1 %>% 
  filter(year >= 1958 & year <= 1962) %>% 
  group_by(standtree) %>% 
  summarise(BAI.pre62 = mean(BAI),
            RD.pre62 = mean(RD),
            LWP.pre62 = mean(LWP),
            EWD.pre62 = mean(EWD),
            LWD.pre62 = mean(LWD),
            MXD.pre62 = mean(MXD)) 

climWide <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", 
                     header = TRUE) %>% 
  select(grow.yr, PDSIavg)

df2 <- df1 %>% 
  mutate(GSL = round(GSL * 0.2296, 1)) %>% 
  mutate(Sample.HT = round(Sample.HT, 1)) %>% 
  left_join(Add58to62, by = "standtree") %>% 
  left_join(climWide, by = c("year" = "grow.yr")) 

n_distinct(df2$standtree)

# Data for climate graphs

PRCP <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", 
                 header = TRUE) %>% 
  select(grow.yr, PRCPtot)

Normals <- data.frame(Month = seq(1,12, by = 1),
                      PRCPmonth = c(2.05, 2.16, 2.12, 1.15, .63, 
                                    .36, 2.61, 3.11, 2.38, 1.66, 
                                    1.76, 1.87)) %>% 
  mutate(PRCPmonth = PRCPmonth * 25.4)

# Plot stats- from TW_Plot_Characteristics.R
plotStats <- read.csv("AnalysisReadyData/TWPaper_plotStats.csv",
                      header = T)

treeStats <- read.csv("AnalysisReadyData/treeData.csv", header = T)

treeSum <- treeStats %>% 
  group_by(PlotID) %>% 
  summarise(DBHsample = mean(DBH.2016),
            SDsample = round(sd(DBH.2016), 2)) %>% 
  mutate(DBHsample = round(DBHsample * 2.54, 2))

plotTreeData <- left_join(plotStats, treeSum, by = "PlotID")
write.csv(plotTreeData, "AnalysisReadyData/plotTreeData.csv", row.names = F)

GSLdata <- plotTreeData %>% 
  group_by(GSL) %>% 
  summarise(BA = mean(BA), DBH = mean(DBH), DBHsd = mean(DBHsd), TPH = mean(TPH), DBHsample = mean(DBHsample),
            SDsample = mean(SDsample))
write.csv(GSLdata, "AnalysisReadyData/GSLdata.csv", row.names = F)

# BA growth plot for Paper

TW_Growth <- read.csv("AnalysisReadyData/TW_Growth_ForPaper.csv", header = T)

# Summarise the pre-1962 stand conditions
TW_sum <- TW_Growth %>% 
  filter(year == 1962 & pre.post == "a") %>% 
  summarise(mean = mean(BA))

df3 <- df2 %>% 
  filter(Sample.HT == 1.4)

# First lines of results
n_distinct(df2$standtree)
n_distinct(df2$TreeID)

mean(df2$RD, na.rm = T)
mean(df2$LWP, na.rm = T)
mean(df2$EWD, na.rm = T)
mean(df2$LWD, na.rm = T)
mean(df2$MXD, na.rm = T)

# Correlations

EWD.RDcor <- df2 %>% 
  na.omit()
cor(EWD.RDcor$EWD, EWD.RDcor$RD) 
cormod <- lm(RD ~ EWD, data = EWD.RDcor)
summary(cormod)

MXD.LWDcor <- df2 %>% 
  na.omit()
cor(MXD.LWDcor$MXD, MXD.LWDcor$LWD)
cor(MXD.LWDcor$RD, MXD.LWDcor$LWP)
cor(MXD.LWDcor$PDSIavg, MXD.LWDcor$EWD)


# Research Question 1- Effect of GSL -------------------------------------

# First, drop samples that don't go back to 1958. Then, filter dataframe to 1963 on. This drops 61 samples, mostly from higher in trees.

df62 <- df2 %>% 
  group_by(standtree) %>% 
  summarise(minyear = min(year)) %>% 
  right_join(df2, by = "standtree") %>% 
  filter(minyear <= 1958) %>% 
  select(-minyear) %>% 
  filter(year >= 1963)

n_distinct(df62$standtree)

# Overall averages- Fig 4

dfsum <- df62 %>% 
  group_by(Sample.HT, GSL, TreeID) %>% 
  summarise(BAImean = mean(BAI, na.rm = T),
            RDmean = mean(RD, na.rm = T),
            LWPmean = mean(LWP, na.rm = T),
            EWDmean = mean(EWD, na.rm = T),
            LWDmean = mean(LWD, na.rm = T),
            MXDmean = mean(MXD, na.rm = T),
            BAIpre = unique(BAI.pre62, na.rm = T),
            RDpre = unique(RD.pre62, na.rm = T),
            LWPpre = unique(LWP.pre62, na.rm = T),
            EWDpre = unique(EWD.pre62, na.rm = T),
            LWDpre = unique(LWD.pre62, na.rm = T),
            MXDpre = unique(MXD.pre62, na.rm = T)) %>% 
  ungroup()

pairwise.fn <- function(x){
  PW1.4 <- CLD(emmeans(x, ~ GSL | Sample.HT, at = list(Sample.HT = 1.4))) %>% 
    select(GSL, PW1.4 = .group)
  PW2.4 <- CLD(emmeans(x, ~ GSL | Sample.HT, at = list(Sample.HT = 2.4))) %>% 
    select(GSL, PW2.4 = .group)
  PW4.9 <- CLD(emmeans(x, ~ GSL | Sample.HT, at = list(Sample.HT = 4.9))) %>% 
    select(GSL, PW4.9 = .group)
  result <- PW1.4 %>% 
    left_join(PW2.4, by = "GSL") %>% 
    left_join(PW4.9, by = "GSL")
  print(result)
}

BAImod <- lme(
  BAImean ~ Sample.HT * factor(GSL) + BAIpre,
  random = ~ 1 | TreeID,
  data = dfsum
)
furnival(BAImod)
Anova(BAImod, type = 3)
pairwise.fn(BAImod)

RDmod <- lme(
  RDmean ~ Sample.HT * factor(GSL) + RDpre,
  random = ~ 1 | TreeID,
  data = dfsum
)
furnival(RDmod)
Anova(RDmod, type = 3)
pairwise.fn(RDmod)

LWPmod <- lme(
  LWPmean ~ Sample.HT * factor(GSL) + LWPpre,
  na.action = na.omit,
  random = ~ 1 | TreeID,
  data = dfsum
)
furnival(LWPmod)
Anova(LWPmod, type = 3)
pairwise.fn(LWPmod)

# Models with year as a covariate

df62 <- df62 %>% 
  na.omit()
  
BAImod <- lme(
  BAI ~ year + PDSIavg + factor(GSL) * Sample.HT + BAI.pre62,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  weights = varPower(form = ~ year),
  data = df62
)
furnival(BAImod)
Anova(BAImod, type = 3)
df62$yhat <- predict(BAImod, level = 0)
bias(df62$yhat, df62$BAI)
pairwise.fn(BAImod)

RDmod <- lme(
  RD ~ year + PDSIavg + factor(GSL) * Sample.HT + RD.pre62,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = df62
)
furnival(RDmod)
Anova(RDmod, type = 3)
df62$yhat <- predict(RDmod, level = 0)
bias(df62$yhat, df62$RD)
pairwise.fn(RDmod)

LWPmod <- lme(
  LWP ~ year + PDSIavg + factor(GSL) * Sample.HT + LWP.pre62,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = df62
)
furnival(LWPmod)
Anova(LWPmod, type = 3)
df62$yhat <- predict(LWPmod, level = 0)
bias(df62$yhat, df62$LWP)
pairwise.fn(LWPmod)

EWDmod <- lme(
  EWD ~ year + PDSIavg + factor(GSL) * Sample.HT + EWD.pre62,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = df62
)
plot(EWDmod)
furnival(EWDmod)
Anova(EWDmod, type = 3)
df62$yhat <- predict(EWDmod, level = 0)
bias(df62$yhat, df62$EWD)

LWDmod <- lme(
  LWD ~ year + PDSIavg + factor(GSL) * Sample.HT + LWD.pre62,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = df62
)
furnival(LWDmod)
Anova(LWDmod, type = 3)
df62$yhat <- predict(LWDmod, level = 0)
bias(df62$yhat, df62$LWD)

MXDmod <- lme(
  MXD ~ year + PDSIavg + factor(GSL) * Sample.HT + MXD.pre62,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = df62
)
furnival(MXDmod)
Anova(MXDmod, type = 3)
df62$yhat <- predict(MXDmod, level = 0)
bias(df62$yhat, df62$MXD)

# Research Question 2- Setting up the data -------------------------------------

# Deal with NAs here
missingRD <- df2 %>%
  group_by(Sample.HT) %>% 
  summarise(Y1951 = mean(RD[year == 1951], na.rm = T),
            Y1956 = mean(RD[year == 1956], na.rm = T),
            Y1996 = mean(RD[year == 1996], na.rm = T),
            Y2002 = mean(RD[year == 2002], na.rm = T)) %>% 
  gather(key = "year", value = "RD.NA", 2:5)

missingMinD <- df2 %>%
  group_by(Sample.HT) %>% 
  summarise(Y1951 = mean(minD[year == 1951], na.rm = T),
            Y1956 = mean(minD[year == 1956], na.rm = T),
            Y1996 = mean(minD[year == 1996], na.rm = T),
            Y2002 = mean(minD[year == 2002], na.rm = T)) %>% 
  gather(key = "year", value = "MinD.NA", 2:5)

missingMXD <- df2 %>%
  group_by(Sample.HT) %>% 
  summarise(Y1951 = mean(MXD[year == 1951], na.rm = T),
            Y1956 = mean(MXD[year == 1956], na.rm = T),
            Y1996 = mean(MXD[year == 1996], na.rm = T),
            Y2002 = mean(MXD[year == 2002], na.rm = T)) %>% 
  gather(key = "year", value = "MXD.NA", 2:5)

missingEWD <- df2 %>%
  group_by(Sample.HT) %>% 
  summarise(Y1951 = mean(EWD[year == 1951], na.rm = T),
            Y1956 = mean(EWD[year == 1956], na.rm = T),
            Y1996 = mean(EWD[year == 1996], na.rm = T),
            Y2002 = mean(EWD[year == 2002], na.rm = T)) %>% 
  gather(key = "year", value = "EWD.NA", 2:5)

missingLWD <- df2 %>%
  group_by(Sample.HT) %>% 
  summarise(Y1951 = mean(LWD[year == 1951], na.rm = T),
            Y1956 = mean(LWD[year == 1956], na.rm = T),
            Y1996 = mean(LWD[year == 1996], na.rm = T),
            Y2002 = mean(LWD[year == 2002], na.rm = T)) %>% 
  gather(key = "year", value = "LWD.NA", 2:5)

missingLWP <- df2 %>%
  group_by(Sample.HT) %>% 
  summarise(Y1951 = mean(LWP[year == 1951], na.rm = T),
            Y1956 = mean(LWP[year == 1956], na.rm = T),
            Y1996 = mean(LWP[year == 1996], na.rm = T),
            Y2002 = mean(LWP[year == 2002], na.rm = T)) %>% 
  gather(key = "year", value = "LWP.NA", 2:5)

missingALL <- missingRD %>% 
  left_join(missingMinD, by = c("Sample.HT", "year")) %>% 
  left_join(missingMXD, by = c("Sample.HT", "year")) %>%
  left_join(missingEWD, by = c("Sample.HT", "year")) %>%
  left_join(missingLWD, by = c("Sample.HT", "year")) %>%
  left_join(missingLWP, by = c("Sample.HT", "year")) %>% 
  separate(year, sep = 1, into = c("drop", "year"), convert = T, remove = T) %>% 
  select(-drop)

df3 <- df2 %>% 
  replace_na(list(RD = -999, minDens = -999, MXD = -999, EWD = -999,
                  LWD = -999, LWP = -999))

df4 <- df3 %>% 
  left_join(missingALL, by = c("Sample.HT", "year"))

df5 <- df4 %>% 
  mutate(RD = ifelse(RD == -999, RD.NA, RD),
         minD = ifelse(minD == -999, MinD.NA, minD),
         MXD = ifelse(MXD == -999, MXD.NA, MXD),
         EWD = ifelse(EWD == -999, EWD.NA, EWD),
         LWD = ifelse(LWD == -999, LWD.NA, LWD),
         LWP = ifelse(LWP == -999, LWP.NA, LWP)) %>% 
  # 6 values ended up as NA in these columns, so I am just replacing with global average
  replace_na(list(EWD = mean(df4$EWD, na.rm = T),
                  LWD = mean(df4$LWD, na.rm = T),
                  LWP = mean(df4$LWP, na.rm = T)))

PRCP <- read.csv("AnalysisReadyData/climLong_TWPaper.csv", header = TRUE) %>% 
  na.omit() %>% 
  select(year, month, PRCP) %>% 
  spread(month, PRCP) %>% 
  rename("JAN" = "1", "FEB" = "2", "MAR" = "3", "APR" = "4", "MAY" = "5",
         "JUN" = "6", "JUL" = "7", "AUG" = "8", "SEP" = "9", "OCT" = "10",
         "NOV" = "11", "DEC" = "12") %>% 
  filter(year > 1917 & year < 2017)

TAVG <- read.csv("AnalysisReadyData/climLong_TWPaper.csv", header = TRUE) %>% 
  na.omit() %>% 
  select(year, month, TAVG) %>% 
  spread(month, TAVG) %>% 
  rename("JAN" = "1", "FEB" = "2", "MAR" = "3", "APR" = "4", "MAY" = "5",
         "JUN" = "6", "JUL" = "7", "AUG" = "8", "SEP" = "9", "OCT" = "10",
         "NOV" = "11", "DEC" = "12") %>% 
  filter(year > 1917 & year < 2017)

# Research Question 2- PRCP at Breast Height --------------------

treeclim.gsl.fn <- function(x,y,gsl1,gsl2){
  
  prechron <- df5 %>%
    filter(Sample.HT == x) %>%
    filter(GSL == gsl1 | GSL == gsl2) %>%
    select(standtree, year, response = y) %>%
    spread(key = standtree, value = response) %>%
    column_to_rownames(var = "year")
  
  prechron.rwi <- detrend(prechron, method = "Spline")
  
  chron <- chron(prechron.rwi, prewhiten = T) %>%
    subset(samp.depth >= 5) %>% 
    rename(index = xxxres) %>% 
    select(index)
  
  climtest <- dcc(
    chrono = chron,
    climate = list(PRCP, TAVG),
    var_names = c("PRCP", "TAVG"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9) +
      .mean("TAVG",-10:-12) +
      .mean("TAVG", 1:3) +
      .mean("TAVG", 4:6) +
      .mean("TAVG", 7:9) +
      .mean("TAVG", -10:9),
    method = "response",
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = x,
           quarter = rep(c("ond", "JFM", "AMJ", "JAS", "Annual"), 2))
}

# RD

BH.RD.low.PRCP <- treeclim.gsl.fn(1.4,"RD",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Low")

BH.RD.mid.PRCP <- treeclim.gsl.fn(1.4,"RD",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
    mutate(group = "Mid")

BH.RD.high.PRCP <- treeclim.gsl.fn(1.4,"RD",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "High")

BH.RD.PRCP <- rbind(BH.RD.low.PRCP, BH.RD.mid.PRCP, BH.RD.high.PRCP) 
BH.RD.PRCP$group <- factor(BH.RD.PRCP$group, levels = c("Low", "Mid", "High"))

# LWP

BH.LWP.low.PRCP <- treeclim.gsl.fn(1.4,"LWP",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Low")

BH.LWP.mid.PRCP <- treeclim.gsl.fn(1.4,"LWP",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Mid")

BH.LWP.high.PRCP <- treeclim.gsl.fn(1.4,"LWP",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "High")

BH.LWP.PRCP <- rbind(BH.LWP.low.PRCP, BH.LWP.mid.PRCP, BH.LWP.high.PRCP) 
BH.LWP.PRCP$group <- factor(BH.LWP.PRCP$group, 
                            levels = c("Low", "Mid", "High"))

# EWD

BH.EWD.low.PRCP <- treeclim.gsl.fn(1.4,"EWD",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Low")

BH.EWD.mid.PRCP <- treeclim.gsl.fn(1.4,"EWD",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Mid")

BH.EWD.high.PRCP <- treeclim.gsl.fn(1.4,"EWD",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "High")

BH.EWD.PRCP <- rbind(BH.EWD.low.PRCP, BH.EWD.mid.PRCP, BH.EWD.high.PRCP) 
BH.EWD.PRCP$group <- factor(BH.EWD.PRCP$group, 
                            levels = c("Low", "Mid", "High"))

# MXD

BH.MXD.low.PRCP <- treeclim.gsl.fn(1.4,"MXD",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Low")

BH.MXD.mid.PRCP <- treeclim.gsl.fn(1.4,"MXD",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "Mid")

BH.MXD.high.PRCP <- treeclim.gsl.fn(1.4,"MXD",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "PRCP") %>% 
  mutate(group = "High")

BH.MXD.PRCP <- rbind(BH.MXD.low.PRCP, BH.MXD.mid.PRCP, BH.MXD.high.PRCP) 
BH.MXD.PRCP$group <- factor(BH.MXD.PRCP$group, 
                            levels = c("Low", "Mid", "High"))

# Final results
BH.RD.PRCP$Response <- "RD"
BH.LWP.PRCP$Response <- "LWP"
BH.EWD.PRCP$Response <- "EWD"
BH.MXD.PRCP$Response <- "MXD"
PRCPresults <- rbind(BH.RD.PRCP, BH.LWP.PRCP, BH.EWD.PRCP, BH.MXD.PRCP)

write.csv(PRCPresults, "RQ2_Results_PRCP.csv", row.names = F)


# Research Question 2- TAVG, at Breast Height --------------------

treeclim.gsl.fn <- function(x,y,gsl1,gsl2){
  
  prechron <- df5 %>%
    filter(Sample.HT == x) %>%
    filter(GSL == gsl1 | GSL == gsl2) %>%
    select(standtree, year, response = y) %>%
    spread(key = standtree, value = response) %>%
    column_to_rownames(var = "year")
  
  prechron.rwi <- detrend(prechron, method = "Spline")
  
  chron <- chron(prechron.rwi, prewhiten = T) %>%
    subset(samp.depth >= 5) %>% 
    rename(index = xxxres) %>% 
    select(index)
  
  climtest <- dcc(
    chrono = chron,
    climate = list(PRCP, TAVG),
    var_names = c("PRCP", "TAVG"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9) +
      .mean("TAVG",-10:-12) +
      .mean("TAVG", 1:3) +
      .mean("TAVG", 4:6) +
      .mean("TAVG", 7:9) +
      .mean("TAVG", -10:9),
    method = "response",
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = x,
           quarter = rep(c("ond", "JFM", "AMJ", "JAS", "Annual"), 2))
}

# RD

BH.RD.low.TAVG <- treeclim.gsl.fn(1.4,"RD",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Low")

BH.RD.mid.TAVG <- treeclim.gsl.fn(1.4,"RD",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Mid")

BH.RD.high.TAVG <- treeclim.gsl.fn(1.4,"RD",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "High")

BH.RD.TAVG <- rbind(BH.RD.low.TAVG, BH.RD.mid.TAVG, BH.RD.high.TAVG) 
BH.RD.TAVG$group <- factor(BH.RD.TAVG$group, levels = c("Low", "Mid", "High"))

# LWP

BH.LWP.low.TAVG <- treeclim.gsl.fn(1.4,"LWP",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Low")

BH.LWP.mid.TAVG <- treeclim.gsl.fn(1.4,"LWP",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Mid")

BH.LWP.high.TAVG <- treeclim.gsl.fn(1.4,"LWP",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "High")

BH.LWP.TAVG <- rbind(BH.LWP.low.TAVG, BH.LWP.mid.TAVG, BH.LWP.high.TAVG) 
BH.LWP.TAVG$group <- factor(BH.LWP.TAVG$group, 
                            levels = c("Low", "Mid", "High"))

# EWD

BH.EWD.low.TAVG <- treeclim.gsl.fn(1.4,"EWD",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Low")

BH.EWD.mid.TAVG <- treeclim.gsl.fn(1.4,"EWD",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Mid")

BH.EWD.high.TAVG <- treeclim.gsl.fn(1.4,"EWD",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "High")

BH.EWD.TAVG <- rbind(BH.EWD.low.TAVG, BH.EWD.mid.TAVG, BH.EWD.high.TAVG) 
BH.EWD.TAVG$group <- factor(BH.EWD.TAVG$group, 
                            levels = c("Low", "Mid", "High"))

# MXD

BH.MXD.low.TAVG <- treeclim.gsl.fn(1.4,"MXD",6.9,13.8) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Low")

BH.MXD.mid.TAVG <- treeclim.gsl.fn(1.4,"MXD",18.4,23) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "Mid")

BH.MXD.high.TAVG <- treeclim.gsl.fn(1.4,"MXD",27.6,34.4) %>% 
  separate(varname, sep = 4, into = c("variable", "drop")) %>% 
  filter(variable == "TAVG") %>% 
  mutate(group = "High")

BH.MXD.TAVG <- rbind(BH.MXD.low.TAVG, BH.MXD.mid.TAVG, BH.MXD.high.TAVG) 
BH.MXD.TAVG$group <- factor(BH.MXD.TAVG$group, 
                            levels = c("Low", "Mid", "High"))

# Final results

BH.RD.TAVG$Response <- "RD"
BH.LWP.TAVG$Response <- "LWP"
BH.EWD.TAVG$Response <- "EWD"
BH.MXD.TAVG$Response <- "MXD"
TAVGresults <- rbind(BH.RD.TAVG, BH.LWP.TAVG, BH.EWD.TAVG, BH.MXD.TAVG)

write.csv(TAVGresults, "RQ2_Results_TAVG.csv", row.names = F)