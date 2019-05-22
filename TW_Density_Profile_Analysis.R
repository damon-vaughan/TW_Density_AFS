library(tidyverse)
library(car)
library(nlme)
library(mgcv)
library(gridExtra)
library(dplR)
library(treeclim)
library(emmeans)
library(dichromat)
source('furnival.r')
source("bias function.r")

# Get dataset ready -------------------------------------------------------

sumProfile <- read.csv("AnalysisReadyData/sumProfileFull_TWPaper_Final.csv", 
                       header = TRUE,
                       stringsAsFactors = FALSE) 

# Drop scans with poor correlations (From crossdating)
dropped <- read.csv("RawData/dropped_scans.csv", header = TRUE, 
                    stringsAsFactors = FALSE) 
dropped <- dropped$Scan

sumProfile2 <- sumProfile[!sumProfile$standtree %in% dropped,]

sumProfile3 <- sumProfile2 %>% 
  filter(standtree != "TW-6-237-45" | year > 1963) %>% 
  filter(standtree != "TW-6-237-1" | year > 1963) %>% 
  filter(standtree != "TW-12-296-45" | year > 1963) %>% 
  filter(standtree != "TW-13-176-1" | (year < 2011 & year > 1948)) %>% 
  filter(standtree != "TW-17-163-1" | year > 1944) 

# Create column to indicate where rings should not be used in intra-ring settings (still ok for ring width and BAI). Criteria are 
# 1) density is less than 200 Results from sample defects such as pitch pockets, decay, and breaks
# 2) LW density is less than EW density. Results from tiny rings that are curved enough that EW does not show up on scan. From B to P, looks like the normal rise for a ring, then a decreased slope, then a radically increased slope where the next ring begins. 

sumProfile4 <- sumProfile3 %>%
  mutate(flawed = ifelse(RD < 200 | is.na(RD), "Y",
                         ifelse(LWD < EWD | is.na(LWD), 
                                "Y", "N")))

# Drop rings with too high a "proportion to min density". Results from compression wood, curved rings, and maybe extractives? Drop them if more than 1x the sd above the mean value.

sumProfile5 <- sumProfile4 %>% 
  mutate(propToMin = distToMin/RW)

cutoff <- mean(sumProfile5$propToMin, na.rm = TRUE) + 
  sd(sumProfile5$propToMin, na.rm = TRUE)

sumProfile6 <- sumProfile5 %>% 
  mutate(flawed = ifelse(propToMin > cutoff | is.na(propToMin), "Y", flawed))

# Latewood proportion and BAI in cm. Area of an annulus is pi(R^2 - r^2). First, adjust mmRingStart by 1 mm to assume half of a 2 mm pith.

sumProfile7 <- sumProfile6 %>% 
  mutate(LWP = LWW/RW,
         mmRingStart = mmRingStart + 1,
         BAI = pi * ((mmRingStart + RW)^2 - mmRingStart^2)/100)

# Filter out first 5 rings and final ring, filter to heights <= 24 feet.

sumProfile8 <- sumProfile7 %>% 
  filter(Sample.HT <= 8) %>%
  group_by(SampleID) %>% 
  filter(newring >= 10 & newring != max(newring)) %>% 
  ungroup() 

# 2 ways to drop based on flawed. Either drop a sample if too many are flawed, or drop individual rings

check <- sumProfile8 %>% 
  group_by(standtree) %>% 
  summarise(total = length(standtree), flawed = length(which(flawed == "Y")), prop = flawed/total)


check <- sumProfile8 %>% 
  group_by(GSL) %>% 
  summarise(count = length(unique(TreeID)))

write.csv(sumProfile8, "AnalysisReadyData/sumProfile20_Final.csv", 
          row.names = FALSE)


# Read in full one --------------------------------------------------------

df1 <- read.csv("AnalysisReadyData/sumProfile20_Final.csv", 
                header = TRUE, stringsAsFactors = FALSE) 

Add58to62 <- df1 %>% 
  filter(year >= 1958 & year <= 1962) %>% 
  group_by(SampleID) %>% 
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
  left_join(Add58to62, by = "SampleID") %>% 
  left_join(climWide, by = c("year" = "grow.yr")) 

n_distinct(df2$SampleID)

# Data Summaries ---------------------------------------------------

# Climate graphs

PRCP <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", 
                 header = TRUE) %>% 
  select(grow.yr, PRCPtot)

Fig1a <- ggplot(PRCP) +
  geom_line(aes(x = grow.yr, y = PRCPtot)) +
  geom_smooth(aes(x = grow.yr, y = PRCPtot), method = "gam", se = FALSE, 
              formula = y ~ s(x, k = 5)) +
  ylab("Total Precipitation (mm)") +
  xlab("Water Year") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))  +
  scale_x_continuous(breaks = seq(1920, 2020, by = 10))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/AnnualPRCP.tiff",
       plot = Fig1a, device = "tiff", 
       width = 25, height = 16, units = "cm")

Fig1a <- ggplot(PRCP) +
  geom_line(aes(x = grow.yr, y = PRCPtot)) +
  geom_smooth(aes(x = grow.yr, y = PRCPtot), method = "gam", se = FALSE, 
              formula = y ~ s(x, k = 5)) +
  ylab("Total Precipitation (mm)") +
  xlab("Water Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))  +
  scale_x_continuous(breaks = seq(1920, 2020, by = 20))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig1a.eps",
       plot = Fig1a, device = "eps", 
       width = 84, height = 60, units = "mm")

# From https://www.ncdc.noaa.gov/cdo-web/datatools/normals it comes in inches

Normals <- data.frame(Month = seq(1,12, by = 1),
                      PRCPmonth = c(2.05, 2.16, 2.12, 1.15, .63, 
                                    .36, 2.61, 3.11, 2.38, 1.66, 
                                    1.76, 1.87)) %>% 
  mutate(PRCPmonth = PRCPmonth * 25.4)

Fig1b <- ggplot(Normals) +
  geom_line(aes(x = Month, y = PRCPmonth)) +
  ylab("Total Precipitation (mm)") +
  xlab("Month") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))  +
  scale_x_continuous(breaks = seq(1,12, by = 1))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/MonthPRCP.tiff",
       plot = Fig1b, device = "tiff", 
       width = 18, height = 11.45, units = "cm")

climWide <- read.csv("AnalysisReadyData/climWide.csv", header = TRUE)

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

Fig2 <- ggplot(TW_Growth) +
  stat_summary(aes(x = year, y = BA, color = factor(GSL)), 
               fun.y = mean, geom = "line", size = 1) +
  geom_hline(yintercept = c(6.9, 13.8, 18.4, 23, 27.6, 34.4), 
             color = colorschemes$Categorical.12[c(1,6,8,10,11,12)],
             linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)) +
  ylab(expression(paste("Basal area ", "(m" ^ "2 ","ha" ^ "-1", ")"))) +
  labs(color = bquote(atop(GSL~phantom(),(m^2~ha^-1)))) +
  scale_x_continuous(breaks = seq(1960, 2020, by = 10)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  guides(color = guide_legend(reverse = TRUE))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig2.tiff",
       device = "tiff", plot = Fig2,
       width = 28, height = 17.5, units = "cm")

# Article figure. Subsetted to BH. Does not look too much different if including all Sample.HT

df3 <- df2 %>% 
  filter(Sample.HT == 1.4)

Fig3a <- ggplot(filter(climWide, grow.yr >= 1940 & grow.yr <= 2016)) +
  geom_line(aes(x = grow.yr, y = PDSIavg)) +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("PDSI") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3b <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = RD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("RD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3c <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = EWD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("EWD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3d <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = LWP), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("LWP") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3e <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = LWD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  ylab("LWD") +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  scale_x_continuous(breaks = seq(1940, 2010, by = 10)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank())

Fig3aGrob <- ggplotGrob(Fig3a)
Fig3bGrob <- ggplotGrob(Fig3b)
Fig3cGrob <- ggplotGrob(Fig3c)
Fig3dGrob <- ggplotGrob(Fig3d)
Fig3eGrob <- ggplotGrob(Fig3e)

Fig3Grob <- rbind(Fig3aGrob, Fig3bGrob, Fig3cGrob, 
                  Fig3dGrob, Fig3eGrob, size = "first")

Fig3 <- grid.arrange(Fig3Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig3.tiff",
       device = "tiff", plot = Fig3,
       width = 28, height = 17.5, units = "cm")

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

Fig4a <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9), 
       aes(x = factor(GSL), y = BAImean, color = factor(Sample.HT))) +
  geom_boxplot(size = .75, notch = F, width = .6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 24, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  ylab(expression(paste("BAI ", "(cm" ^ "2", ")"))) +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

Fig4b <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9), 
       aes(x = factor(GSL), y = RDmean, color = factor(Sample.HT))) +
  geom_boxplot(size = .75, notch = F, width = .6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 24, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  ylab(expression(paste("RD ", "(kg m" ^ "-3", ")"))) +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

Fig4c <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9),  
       aes(x = factor(GSL), y = LWPmean, color = factor(Sample.HT))) +
  geom_boxplot(size = .75, notch = F, width = .6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 24, margin = margin(r = 15)),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  xlab(expression(paste("GSL ", "(m" ^ "2 ","ha" ^ "-1", ")"))) +
  ylab("LWP") +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

Fig4aGrob <- ggplotGrob(Fig4a)
Fig4bGrob <- ggplotGrob(Fig4b)
Fig4cGrob <- ggplotGrob(Fig4c)
Fig4Grob <- rbind(Fig4aGrob, Fig4bGrob, Fig4cGrob, size = "first")
Fig4 <- grid.arrange(Fig4Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig4.tiff", plot = Fig4,
       device = "tiff", 
       width = 28, height = 17.5, units = "cm")

Fig4legend <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9),   
                     aes(x = factor(GSL), y = BAImean)) +
  geom_boxplot(aes(color = factor(Sample.HT)), size = .75) +
  theme_bw() +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)) +
  labs(color = "Sample\nheight (m)") +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig4Legend.tiff", 
       plot = Fig4legend,
       device = "tiff", 
       width = 28, height = 17.5, units = "cm")

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
plot(BAImod)
furnival(BAImod)
Anova(BAImod, type = 3)
# summary(BAImod)$coefficients$fixed
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
plot(RDmod)
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
plot(LWPmod)
furnival(LWPmod)
Anova(LWPmod, type = 3)
# summary(LWPmod)$coefficients$fixed
df62$yhat <- predict(LWPmod, level = 0)
bias(df62$yhat, df62$LWP)
pairwise.fn(LWPmod)
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

Fig5a <- ggplot(filter(df62, Sample.HT == 1.4), 
       aes(x = year, y = BAI, color = factor(GSL))) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), se = F, size = 1) +
  ylab(expression(paste("BAI ", "(cm" ^ "2", ")"))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 24, margin = margin(r = 15))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(1963, 2016), breaks = seq(1970, 2010, by = 10)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  theme(plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "cm"))

Fig5b <- ggplot(filter(df62, Sample.HT == 1.4), 
       aes(x = year, y = RD, color = factor(GSL))) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), se = F, size = 1) +
  ylab(expression(paste("RD ", "(kg m" ^ "-3", ")"))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 24, margin = margin(r = 15))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(1963, 2016), breaks = seq(1970, 2010, by = 10)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  theme(plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "cm"))

Fig5c <- ggplot(filter(df62, Sample.HT == 1.4),  
       aes(x = year, y = LWP, color = factor(GSL))) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), se = F, size = 1) +
  ylab("LWP") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 24, margin = margin(r = 15))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(1963, 2016), breaks = seq(1970, 2010, by = 10)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  theme(plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "cm"))

Fig5aGrob <- ggplotGrob(Fig5a)
Fig5bGrob <- ggplotGrob(Fig5b)
Fig5cGrob <- ggplotGrob(Fig5c)
Fig5Grob <- rbind(Fig5aGrob, Fig5bGrob, Fig5cGrob, size = "first")
Fig5 <- grid.arrange(Fig5Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig5.tiff", plot = Fig5,
       device = "tiff", 
       width = 16, height = 20, units = "cm")

Fig5Legend <- ggplot(df62) +
  geom_line(aes(x = year, y = BAI, color = factor(GSL)), size = 1) +
  theme_bw() +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)) +
  labs(color = bquote(atop(GSL~phantom(),(m^2~ha^-1)))) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)])

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig5_Legend.tiff", plot = Fig5Legend,
       device = "tiff", 
       width = 28, height = 17.5, units = "cm")



# Research Question 2- Getting set up -------------------------------------

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
  # 6 random values ended up as NA in these columns, so I am just replacing with global average
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

# Research Question 2- GSL separation, at BH, for PRCP --------------------


treeclim.gsl.fn <- function(x,y,gsl1,gsl2){
  
  prechron <- df5 %>%
    filter(Sample.HT == x) %>%
    filter(GSL == gsl1 | GSL == gsl2) %>%
    select(SampleID, year, response = y) %>%
    spread(key = SampleID, value = response) %>%
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

Fig6a <- ggplot(BH.RD.PRCP,
       aes(x = quarter, y = coef, ymin = ci_lower, 
           ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.RD.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  ylab("Correlation coefficient") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  ggtitle("RD") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 24)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

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

Fig6b <- ggplot(BH.LWP.PRCP,
       aes(x = quarter, y = coef, ymin = ci_lower, 
           ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.LWP.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  ggtitle("LWP") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 24)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

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

Fig6c <- ggplot(BH.EWD.PRCP,
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.EWD.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  ggtitle("EWD") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 24)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

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

Fig6d <- ggplot(BH.MXD.PRCP,
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.MXD.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  ggtitle("MXD") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 24)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

Fig6aGrob <- ggplotGrob(Fig6a)
Fig6bGrob <- ggplotGrob(Fig6b)
Fig6cGrob <- ggplotGrob(Fig6c)
Fig6dGrob <- ggplotGrob(Fig6d)

Fig6_1Grob <- cbind(Fig6aGrob, Fig6bGrob, Fig6cGrob, Fig6dGrob, size = "first")

Fig6.1 <- grid.arrange(Fig6_1Grob)

Fig6.1Legend <- ggplot(BH.EWD.PRCP,
                     aes(x = quarter, y = coef, ymin = ci_lower, 
                         ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  ylab("Precipitation\nResponse coef.") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20)) +
  labs(color = "GSL group") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)) +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig6_1Legend.tiff", plot = Fig6.1Legend,
       device = "tiff", 
       width = 28, height = 7.5, units = "cm")

# Final results
BH.RD.PRCP$Response <- "RD"
BH.LWP.PRCP$Response <- "LWP"
BH.EWD.PRCP$Response <- "EWD"
BH.MXD.PRCP$Response <- "MXD"
PRCPresults <- rbind(BH.RD.PRCP, BH.LWP.PRCP, BH.EWD.PRCP, BH.MXD.PRCP)

write.csv(PRCPresults, "RQ2_Results_PRCP.csv", row.names = F)



# Research Question 2- GSL separation, at BH, for TAVG --------------------

treeclim.gsl.fn <- function(x,y,gsl1,gsl2){
  
  prechron <- df5 %>%
    filter(Sample.HT == x) %>%
    filter(GSL == gsl1 | GSL == gsl2) %>%
    select(SampleID, year, response = y) %>%
    spread(key = SampleID, value = response) %>%
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

Fig7a <- ggplot(BH.RD.TAVG,
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.RD.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

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

Fig7b <- ggplot(BH.LWP.TAVG,
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.LWP.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

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

Fig7c <- ggplot(BH.EWD.TAVG,
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.EWD.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

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

Fig7d <- ggplot(BH.MXD.TAVG,
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete(limits = BH.MXD.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(size = 20)) +
  theme(legend.position = "none") +
  scale_color_manual(values = colorschemes$Categorical.12[c(1,10,12)])

Fig7aGrob <- ggplotGrob(Fig7a)
Fig7bGrob <- ggplotGrob(Fig7b)
Fig7cGrob <- ggplotGrob(Fig7c)
Fig7dGrob <- ggplotGrob(Fig7d)

Fig6_2Grob <- cbind(Fig7aGrob, Fig7bGrob, Fig7cGrob, Fig7dGrob, size = "first")

Fig6Grob <- rbind(Fig6_1Grob, Fig6_2Grob, size = "first")

Fig6 <- grid.arrange(Fig6Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig6.tiff", plot = Fig6,
       device = "tiff", 
       width = 28, height = 17.5, units = "cm")

# legend

Fig6_2Legend <- ggplot(BH.EWD.TAVG,
       aes(x = quarter, y = coef, ymin = ci_lower, 
           ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = 1) +
  ylab("Temperature\nResponse coef.") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 20)) +
  labs(color = "GSL group") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig6_2Legend.tiff", plot = Fig6_2Legend,
       device = "tiff", 
       width = 28, height = 17.5, units = "cm")

# Final results

BH.RD.TAVG$Response <- "RD"
BH.LWP.TAVG$Response <- "LWP"
BH.EWD.TAVG$Response <- "EWD"
BH.MXD.TAVG$Response <- "MXD"
TAVGresults <- rbind(BH.RD.TAVG, BH.LWP.TAVG, BH.EWD.TAVG, BH.MXD.TAVG)

write.csv(TAVGresults, "RQ2_Results_TAVG.csv", row.names = F)

# Research Question 2- No GSL separation ----------------------------------

treeclim.fn <- function(x,y){
  
  prechron <- df5 %>%
    filter(Sample.HT == x) %>% 
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
    climate = list(PRCP, TAVG, TMIN, TMAX),
    var_names = c("PRCP", "TAVG", "TMIN", "TMAX"),
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

climtest <- treeclim.fn(1.4,9)

ringDensplot <- ggplot(climtest, 
                       aes(x = quarter, y = coef, ymin = ci_lower, 
                           ymax = ci_upper, color = factor(varname))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = climtest$quarter[1:5]) +
  ggtitle("Ring density") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Climvar")

# Results at BH

BAIResults <- treeclim.fn(1.4, 32) %>% 
  mutate(response = "BAI")

BAIplot <- ggplot(filter(BAIResults, varname == "PRCP.sum"), 
                  aes(x = quarter, y = coef, ymin = ci_lower, 
                      ymax = ci_upper, color = factor(varname))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = BAIResults$quarter[1:5]) +
  ggtitle("BAI") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Climvar")

ringDensResults <- treeclim.fn(1.4, 9) %>% 
  mutate(response = "ringDens")

ringDensplot <- ggplot(filter(ringDensResults, varname == "PRCP.sum"), 
                       aes(x = quarter, y = coef, ymin = ci_lower, 
                           ymax = ci_upper, color = factor(varname))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = ringDensResults$quarter[1:5]) +
  ggtitle("RD") +
  ylab("Correlation") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Climvar") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = F)

LWpropResults <- treeclim.fn(1.4, 31) %>% 
  mutate(response = "LWprop")

LWpropplot <- ggplot(filter(LWpropResults, varname == "PRCP.sum"), 
                     aes(x = quarter, y = coef, ymin = ci_lower, 
                         ymax = ci_upper, color = factor(varname))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = LWpropResults$quarter[1:5]) +
  ggtitle("LWP") +
  ylab("Correlation") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Climvar") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = F)

EWdensResults <- treeclim.fn(1.4, 12) %>% 
  mutate(response = "EWdens")

EWdensplot <- ggplot(filter(EWdensResults, varname == "PRCP.sum"), 
                     aes(x = quarter, y = coef, ymin = ci_lower, 
                         ymax = ci_upper, color = factor(varname))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = EWdensResults$quarter[1:5]) +
  ggtitle("EWD") +
  ylab("Correlation") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Climvar") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = F)

maxDensResults <- treeclim.fn(1.4, 11) %>% 
  mutate(response = "maxDens")

maxDensplot <- ggplot(filter(maxDensResults, varname == "PRCP.sum"), 
                      aes(x = quarter, y = coef, ymin = ci_lower, 
                          ymax = ci_upper, color = factor(varname))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = maxDensResults$quarter[1:5]) +
  ggtitle("MXD") +
  ylab("Correlation") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Climvar") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = F)

ggarrange(ringDensplot, LWpropplot, EWdensplot, maxDensplot) 

# The other way to present- each climate variable on separate plot, compare responses next to each other

Results <- rbind(ringDensResults, LWpropResults, 
                 EWdensResults, maxDensResults)

# Plot of all responses for PRCP
ggplot(filter(Results, varname == "PRCP.sum"), 
       aes(x = quarter, y = coef, ymin = ci_lower, 
           ymax = ci_upper, color = factor(response))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = LWdensResults$quarter[1:5]) +
  ggtitle("PRCP") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Density component")

# Plot of all responses for TAVG
ggplot(filter(Results, varname == "TAVG.mean"), 
       aes(x = quarter, y = coef, ymin = ci_lower, 
           ymax = ci_upper, color = factor(response))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = LWdensResults$quarter[1:5]) +
  ggtitle("TAVG") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Density component")  








##################################
# Old Stuff!!!!


# Research Question 3- with chronologies ------------------------------------

climWide <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", header = TRUE)

YSTstart <- data.frame(year = seq(1918, 1961), YST = rep (0, 44))
YSTmid <- data.frame(year = seq(1962, 2001), YST = rep(seq(0,9), 4))
YSTmid2 <- data.frame(year = 2002, YST = 10)
YSTend <- data.frame(year = seq(2003, 2016), YST = seq(0,13))

YSTfull <- rbind(YSTstart,YSTmid,YSTmid2,YSTend)

x <- 4.5
y <- 9

chron.fn <- function(x,y){
  
  missYearAvg <- df2 %>%
    filter(Sample.HT == x) %>%
    group_by(year) %>%
    select(year, response = y) %>% 
    summarise(mean = mean(response, na.rm = T)) %>%
    filter(year == 1951 | year == 1956 | year == 1996 | year == 2002)
  
  prechron <- df2 %>%
    filter(Sample.HT == x) %>% 
    select(standtree, year, response = y) %>%
    replace_na(list(response = mean(missYearAvg$mean))) %>%
    spread(key = standtree, value = response) %>%
    column_to_rownames(var = "year")
  
  prechron.rwi <- detrend(prechron, method = "Spline")
  
  chron <- chron(prechron.rwi, prewhiten = T) %>%
    rename(index = xxxres) %>% 
    select(index) 
}

RDchron <- chron.fn(4.5, 9)
MXDchron <- chron.fn(4.5, 11)
LWPchron <- chron.fn(4.5, 31)
str(df2)
chron2 <- RDchron %>% 
  rownames_to_column(var = "year")
chron2$year <- as.numeric(chron2$year)

chron3 <- chron2 %>% 
  left_join(climWide, by = c("year" = "grow.yr")) %>% 
  left_join(YSTfull) %>% 
  na.omit()

modtest <- lm(index ~ PRCP.2 + TAVG.3 + TAVG.4, data = chron3)
modtest <- gam(index ~ s(PRCP.2) + s(TAVG.3) + s(TAVG.4), data = chron3)
modtest <- gam(index ~ s(PRCP.1) + YST, data = chron3)

modtest <- gam(index ~ s(PDSIavg) + YST, data = chron3)

chron3$yhat = predict(modtest)
summary(modtest)
Anova(modtest, type = 3)

ggplot(chron3) +
  geom_line(aes(x = year, y = yhat)) +
  geom_point(aes(x = year, y = index)) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003))

chron3$resid <- modtest$residuals

ggplot(chron3) +
  geom_line(aes(x = year, y = resid)) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003))

chron362pre <- filter(chron3, year >= 1957 & year < 1962)
chron362post <- filter(chron3, year > 1963 & year <= 1968)
t.test(chron362pre$resid, chron362post$resid)



# Old- Research Question 3- 1962 thinning -----------------------------

library(rgl)

df3 <- df2

climWide <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", 
                     header = TRUE) %>% 
  select(grow.yr, PDSIavg) %>% 
  mutate(dummyrep = rnorm(length(grow.yr), -.6,2.27))

results.fn <- function(x){
  significance <- Anova(x, type = 3) %>% 
    select(-Chisq, -Df) %>% 
    rownames_to_column() %>% 
    rename("covariate" = rowname, "pval" = "Pr(>Chisq)") %>% 
    mutate(pval = ifelse(pval < .001, "***",
                         ifelse(pval < .01, "**",
                                ifelse(pval < .05, "*",
                                       ifelse(pval < .1, ".", " ")))))
  results <- summary(x)$coefficients$fixed %>% 
    enframe() %>% 
    rename("covariate" = name, "effect" = value) %>% 
    mutate(effect = round(effect, 4)) %>% 
    right_join(significance, by = "covariate") %>% 
    gather("type", "value", 2:3) %>% 
    arrange(covariate)
  return(results)
}

YSTstart <- data.frame(year = seq(1918, 1961), YST = rep (0, 44))
YSTmid <- data.frame(year = seq(1962, 2001), YST = rep(seq(0,9), 4))
YSTmid2 <- data.frame(year = 2002, YST = 10)
YSTend <- data.frame(year = seq(2003, 2016), YST = seq(0,13))

YSTfull <- rbind(YSTstart,YSTmid,YSTmid2,YSTend)

# 3d plots

# All years all data
df4 <- df3 %>% 
  left_join(YSTfull, by = "year") %>% 
  left_join(climWide, by = c("year" = "grow.yr")) %>% 
  select(ringDens, YST, PDSIavg, year)
plot3d(df4[,1:3])

# All years means
df4mean <- df4 %>% 
  group_by(YST, PDSIavg) %>% 
  summarise(ringDens = mean(ringDens, na.rm = T))
plot3d(df4mean[,1:3])

# First thinning all data
dfsub <- df4 %>%
  filter(year >= 1958 & year <= 1967)
plot3d(dfsub[,1:3])

# First thinning means
dfsubmean <- dfsub %>% 
  group_by(YST, PDSIavg) %>% 
  summarise(ringDens = mean(ringDens, na.rm = T))
plot3d(dfsubmean[,1:3])

# The 1962 thinning

df4 <- df3 %>% 
  left_join(YSTfull, by = "year") %>% 
  left_join(climWide, by = c("year" = "grow.yr")) %>% 
  mutate(dummyfull = rnorm(length(standtree), -.6,2.27))

df4$Sample.HT <- as.factor(df4$Sample.HT)

Add58to62 <- df4 %>% 
  filter(year >= 1958 & year <= 1962) %>% 
  group_by(standtree) %>% 
  summarise(BAI.pre62 = mean(BAI, na.rm = TRUE),
            ringDens.pre62 = mean(ringDens, na.rm = TRUE),
            LW.prop.pre62 = mean(LW.prop, na.rm = TRUE),
            LW.dens.pre62 = mean(LW.dens, na.rm = TRUE),
            EW.dens.pre62 = mean(EW.dens, na.rm = TRUE)) 

fullyear <- df4 %>%
  left_join(Add58to62, by = "standtree") %>%
  na.omit()

# post10 <- df4 %>% 
#   left_join(Add58to62, by = "standtree") %>% 
#   filter(year >= 1962 & year <= 1972) %>% 
#   na.omit()
# 
# post5 <- df4 %>% 
#   left_join(Add58to62, by = "standtree") %>% 
#   filter(year >= 1962 & year <= 1967) %>% 
#   na.omit()

# NEW- full years
# Three variables that repeat by year- YST, dummyrep, and PDSIavg
# Of the three, by itself, dummyrep is not significant. Add YST, still not. Add PDSIavg, and dummyrep becomes significant. Include all 3, and all three are significant

fullyearmod <- lme(
  ringDens ~
    PDSIavg,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = fullyear
)
Anova(fullyearmod, type = 3)
summary(fullyearmod)
ggplot(fullyear) +
  stat_summary(aes(x = PDSIavg, y = ringDens), fun.y = mean, geom = "point")

plot(fullyear$dummyrep,fullyear$ringDens)

fullyearmod <- lm(
  ringDens ~ Sample.HT * ringDens.pre62 + 
    GSL + dummyrep + dummyfull,
  # random = ~ 1 | TreeID,
  na.action = na.omit,
  data = fullyear
)
Anova(fullyearmod, type = 3)

# Below is same but filtered to post-1962. Any combination explored above results in significant dummyrep
fullyearmod <- lme(
  ringDens ~ Sample.HT * ringDens.pre62 + 
    GSL + dummyrep + PDSIavg + YST + dummyfull,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = filter(fullyear, year >= 1962)
)
Anova(fullyearmod, type = 3)




# BAI response
RQ1.1a <- lme(
  BAI ~ Sample.HT * YST + PDSIavg + BAI.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.1b <- lme(
  BAI ~ Sample.HT * YST + PDSIavg + BAI.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.1c <- lme(
  BAI ~ YST + Sample.HT + PDSIavg + BAI.pre62 + GSL + dummyrep,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
RQ1.1c <- lme(
  dummyfull ~ YST + Sample.HT + PDSIavg + BAI.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)

post10$yhat <- predict(RQ1.1a, level = 0)
bias(post10$yhat, post10$BAI)

post5$yhat <- predict(RQ1.1b, level = 0)
bias(post5$yhat, post5$BAI)

post5$yhat <- predict(RQ1.1c, level = 0)
bias(post5$yhat, post5$BAI)

furnival(RQ1.1c)
summary(RQ1.1c)
Anova(RQ1.1c, type = 3)

results1.1a <- results.fn(RQ1.1a) %>% 
  rename("OG" = value)
results1.1b <- results.fn(RQ1.1b) %>% 
  rename("OGpost5" = value)
results1.1c <- results.fn(RQ1.1c) %>% 
  rename("Post5CAR" = value)
summary(RQ1.1c)
plot(RQ1.1c)

vec <- c("NA","NA", furnival(RQ1.1a)$coefd[2,1],
         furnival(RQ1.1b)$coefd[2,1],
         furnival(RQ1.1c)$coefd[2,1])

results1.1 <- left_join(results1.1a, results1.1b) %>% 
  left_join(results1.1c) %>%
  rbind(vec) 

RQ1.1d <- RQ1.1c
summary(RQ1.1d)
furnival(RQ1.1d)
post5$yhat <- predict(RQ1.1d, level = 0)
bias(post5$yhat, post5$BAI)

ggplot(df4) +
  # geom_smooth(aes(x=year, y = yhat.f, color = factor(Sample.HT))) +
  geom_smooth(aes(x = year, y = yhat.r, color = factor(Sample.HT)))

# Ring density response
RQ1.2a <- lme(
  ringDens ~ Sample.HT * YST + PDSIavg + ringDens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.2b <- lme(
  ringDens ~ Sample.HT * YST + PDSIavg + ringDens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.2c <- lme(
  ringDens ~ Sample.HT * YST + PDSIavg + ringDens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.2c)

post10$yhat <- predict(RQ1.2a, level = 0)
bias(post10$yhat, post10$ringDens)

post5$yhat <- predict(RQ1.2b, level = 0)
bias(post5$yhat, post5$ringDens)

post5$yhat <- predict(RQ1.2c, level = 0)
bias(post5$yhat, post5$ringDens)

results1.2a <- results.fn(RQ1.2a) %>% 
  rename("OG" = value)
results1.2b <- results.fn(RQ1.2b) %>% 
  rename("OGpost5" = value)
results1.2c <- results.fn(RQ1.2c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.2a)$coefd[2,1],
         furnival(RQ1.2b)$coefd[2,1],
         furnival(RQ1.2c)$coefd[2,1])

results1.2 <- left_join(results1.2a, results1.2b) %>% 
  left_join(results1.2c) %>% 
  rbind(vec) 

RQ1.2d <- update(RQ1.2c, ringDens ~ Sample.HT + YST + PDSIavg + ringDens.pre62)
summary(RQ1.2d)
Anova(RQ1.2d, type = 3)
furnival(RQ1.2d)
post5$yhat <- predict(RQ1.2d, level = 0)
bias(post5$yhat, post5$ringDens)

# Latewood proportion response
RQ1.3a <- lme(
  LW.prop ~ Sample.HT * YST + PDSIavg + LW.prop.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.3b <- lme(
  LW.prop ~ Sample.HT * YST + PDSIavg + LW.prop.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.3c <- lme(
  LW.prop ~ Sample.HT * YST + PDSIavg + LW.prop.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.3c)

post10$yhat <- predict(RQ1.3a, level = 0)
bias(post10$yhat, post10$LW.prop)

post5$yhat <- predict(RQ1.3b, level = 0)
bias(post5$yhat, post5$LW.prop)

post5$yhat <- predict(RQ1.3c, level = 0)
bias(post5$yhat, post5$LW.prop)

results1.3a <- results.fn(RQ1.3a) %>% 
  rename("OG" = value)
results1.3b <- results.fn(RQ1.3b) %>% 
  rename("OGPost5" = value)
results1.3c <- results.fn(RQ1.3c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.3a)$coefd[2,1],
         furnival(RQ1.3b)$coefd[2,1],
         furnival(RQ1.3c)$coefd[2,1])

results1.3 <- left_join(results1.3a, results1.3b) %>% 
  left_join(results1.3c) %>% 
  rbind(vec) 

RQ1.3d <- update(RQ1.3c, LW.prop ~ Sample.HT * YST + PDSIavg + LW.prop.pre62)  
Anova(RQ1.3d, type = 3)
summary(RQ1.3d)
furnival(RQ1.3d)

post5$yhat <- predict(RQ1.3d, level = 0)
bias(post5$yhat, post5$LW.prop)

# Earlywood response
RQ1.4a <- lme(
  EW.dens ~ Sample.HT * YST + PDSIavg + EW.dens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.4b <- lme(
  EW.dens ~ Sample.HT * YST + PDSIavg + EW.dens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.4c <- lme(
  EW.dens ~ Sample.HT * YST + PDSIavg + EW.dens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.4c)

post10$yhat <- predict(RQ1.4a, level = 0)
bias(post10$yhat, post10$EW.dens)

post5$yhat <- predict(RQ1.4b, level = 0)
bias(post5$yhat, post5$EW.dens)

post5$yhat <- predict(RQ1.4c, level = 0)
bias(post5$yhat, post5$EW.dens)

results1.4a <- results.fn(RQ1.4a) %>% 
  rename("OG" = value)
results1.4b <- results.fn(RQ1.4b) %>% 
  rename("OGPost5" = value)
results1.4c <- results.fn(RQ1.4c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.4a)$coefd[2,1],
         furnival(RQ1.4b)$coefd[2,1],
         furnival(RQ1.4c)$coefd[2,1])

results1.4 <- left_join(results1.4a, results1.4b) %>% 
  left_join(results1.4c) %>% 
  rbind(vec) 

RQ1.4d <- update(RQ1.4c,  EW.dens ~ Sample.HT * YST + PDSIavg + EW.dens.pre62)
Anova(RQ1.4d, type = 3)
summary(RQ1.4d)
furnival(RQ1.4d)
post5$yhat <- predict(RQ1.4d, level = 0)
bias(post5$yhat, post5$EW.dens)

# Latewood density 
RQ1.5a <- lme(
  LW.dens ~ Sample.HT * YST + PDSIavg + LW.dens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.5b <- lme(
  LW.dens ~ Sample.HT * YST + PDSIavg + LW.dens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.5c <- lme(
  LW.dens ~ Sample.HT * YST + PDSIavg + LW.dens.pre62 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.5c)

post10$yhat <- predict(RQ1.5a, level = 0)
bias(post10$yhat, post10$LW.dens)

post5$yhat <- predict(RQ1.5b, level = 0)
bias(post5$yhat, post5$LW.dens)

post5$yhat <- predict(RQ1.5c, level = 0)
bias(post5$yhat, post5$LW.dens)

results1.5a <- results.fn(RQ1.5a) %>% 
  rename("OG" = value)
results1.5b <- results.fn(RQ1.5b) %>% 
  rename("OGPost5" = value)
results1.5c <- results.fn(RQ1.5c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.5a)$coefd[2,1],
         furnival(RQ1.5b)$coefd[2,1],
         furnival(RQ1.5c)$coefd[2,1])

results1.5 <- left_join(results1.5a, results1.5b) %>% 
  left_join(results1.5c) %>% 
  rbind(vec) 

RQ1.5d <- update(RQ1.5c,  LW.dens ~ YST + PDSIavg + LW.dens.pre62)
Anova(RQ1.5d, type = 3)
summary(RQ1.5d)
furnival(RQ1.5d)
post5$yhat <- predict(RQ1.5d, level = 0)
bias(post5$yhat, post5$LW.dens)

# Old Research question 3- Light thinnings ------------------------------------
results.fn <- function(x){
  significance <- Anova(x, type = 3) %>% 
    select(-Chisq, -Df) %>% 
    rownames_to_column() %>% 
    rename("covariate" = rowname, "pval" = "Pr(>Chisq)") %>% 
    mutate(pval = ifelse(pval < .001, "***",
                         ifelse(pval < .01, "**",
                                ifelse(pval < .05, "*",
                                       ifelse(pval < .1, ".", " ")))))
  results <- summary(x)$coefficients$fixed %>% 
    enframe() %>% 
    rename("covariate" = name, "effect" = value) %>% 
    mutate(effect = round(effect, 4)) %>% 
    right_join(significance, by = "covariate") %>% 
    gather("type", "value", 2:3) %>% 
    arrange(covariate)
  return(results)
}

AddPre1972 <- df3 %>% 
  filter(year >= 1968 & year <= 1972) %>% 
  group_by(standtree) %>% 
  summarise(BAI.pre72 = mean(BAI, na.rm = TRUE),
            ringDens.pre72 = mean(ringDens, na.rm = TRUE),
            LW.prop.pre72 = mean(LW.prop, na.rm = TRUE),
            LW.dens.pre72 = mean(LW.dens, na.rm = TRUE),
            EW.dens.pre72 = mean(EW.dens, na.rm = TRUE)) 

post10 <- df3 %>% 
  left_join(AddPre1972, by = "standtree") %>% 
  filter(year > 1972 & YST <= 9) %>% 
  na.omit()

post5 <- df3 %>% 
  left_join(AddPre1972, by = "standtree") %>% 
  filter(year > 1972 & YST <= 5) %>% 
  na.omit()

# BAI Response after 1972
RQ1.6a <- lme(
  BAI ~ YST * Sample.HT + PDSIavg + BAI.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.6b <- lme(
  BAI ~ YST * Sample.HT + PDSIavg + BAI.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.6c <- lme(
  BAI ~ YST * Sample.HT + PDSIavg + BAI.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.6c)

post10$yhat <- predict(RQ1.6a, level = 0)
bias(post10$yhat, post10$BAI)

post5$yhat <- predict(RQ1.6b, level = 0)
bias(post5$yhat, post5$BAI)

post5$yhat <- predict(RQ1.6c, level = 0)
bias(post5$yhat, post5$BAI)

results1.6a <- results.fn(RQ1.6a) %>% 
  rename("OG" = value)
results1.6b <- results.fn(RQ1.6b) %>% 
  rename("OGPost5" = value)
results1.6c <- results.fn(RQ1.6c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.6a)$coefd[2,1],
         furnival(RQ1.6b)$coefd[2,1],
         furnival(RQ1.6c)$coefd[2,1])

results1.6 <- left_join(results1.6a, results1.6b) %>% 
  left_join(results1.6c) %>%
  rbind(vec) 

RQ1.6d <- update(RQ1.6c, BAI ~ YST + PDSIavg + BAI.pre72 + GSL)
Anova(RQ1.6d)
summary(RQ1.6d)
furnival(RQ1.6d)
post5$yhat <- predict(RQ1.6d, level = 0)
bias(post5$yhat, post5$BAI)

# Ring density response post 1972
RQ1.7a <- lme(
  ringDens ~ YST * Sample.HT + PDSIavg + ringDens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.7b <- lme(
  ringDens ~ YST * Sample.HT + PDSIavg + ringDens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.7c <- lme(
  ringDens ~ YST * Sample.HT + PDSIavg + ringDens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.7c)

results1.7a <- results.fn(RQ1.7a) %>% 
  rename("OG" = value)
results1.7b <- results.fn(RQ1.7b) %>% 
  rename("OGPost5" = value)
results1.7c <- results.fn(RQ1.7c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.7a)$coefd[2,1],
         furnival(RQ1.7b)$coefd[2,1],
         furnival(RQ1.7c)$coefd[2,1])

results1.7 <- left_join(results1.7a, results1.7b) %>% 
  left_join(results1.7c) %>%
  rbind(vec) 

RQ1.7d <- update(RQ1.7c, ringDens ~ YST * Sample.HT + PDSIavg + ringDens.pre72)
Anova(RQ1.7d, type = 3)
summary(RQ1.7d)
furnival(RQ1.7d)
post5$yhat <- predict(RQ1.7d, level = 0)
bias(post5$yhat, post5$ringDens)

# Latewood proportion post 1972
RQ1.8a <- lme(
  LW.prop ~ YST * Sample.HT + PDSIavg + LW.prop.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.8b <- lme(
  LW.prop ~ YST * Sample.HT + PDSIavg + LW.prop.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.8c <- lme(
  LW.prop ~ YST * Sample.HT + PDSIavg + LW.prop.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.8c)

post10$yhat <- predict(RQ1.8a, level = 0)
bias(post10$yhat, post10$LW.prop)

post5$yhat <- predict(RQ1.8b, level = 0)
bias(post5$yhat, post5$LW.prop)

post5$yhat <- predict(RQ1.8c, level = 0)
bias(post5$yhat, post5$LW.prop)

results1.8a <- results.fn(RQ1.8a) %>% 
  rename("OG" = value)
results1.8b <- results.fn(RQ1.8b) %>% 
  rename("OGPost5" = value)
results1.8c <- results.fn(RQ1.8c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.8a)$coefd[2,1],
         furnival(RQ1.8b)$coefd[2,1],
         furnival(RQ1.8c)$coefd[2,1])

results1.8 <- left_join(results1.8a, results1.8b) %>% 
  left_join(results1.8c) %>%
  rbind(vec) 

RQ1.8d <- update(RQ1.8c, LW.prop ~ YST + Sample.HT + PDSIavg + LW.prop.pre72 + GSL)
Anova(RQ1.8d)
summary(RQ1.8d)
furnival(RQ1.8d)
post5$yhat <- predict(RQ1.8d, level = 0)
bias(post5$yhat, post5$LW.prop)

# EW dens post 1972
RQ1.9a <- lme(
  EW.dens ~ YST * Sample.HT + PDSIavg + EW.dens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.9b <- lme(
  EW.dens ~ YST * Sample.HT + PDSIavg + EW.dens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.9c <- lme(
  EW.dens ~ YST * Sample.HT + PDSIavg + EW.dens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.9c)

post10$yhat <- predict(RQ1.9a, level = 0)
bias(post10$yhat, post10$EW.dens)

post5$yhat <- predict(RQ1.9b, level = 0)
bias(post5$yhat, post5$EW.dens)

post5$yhat <- predict(RQ1.9c, level = 0)
bias(post5$yhat, post5$EW.dens)

results1.9a <- results.fn(RQ1.9a) %>% 
  rename("OG" = value)
results1.9b <- results.fn(RQ1.9b) %>% 
  rename("OGPost5" = value)
results1.9c <- results.fn(RQ1.9c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.9a)$coefd[2,1],
         furnival(RQ1.9b)$coefd[2,1],
         furnival(RQ1.9c)$coefd[2,1])

results1.9 <- left_join(results1.9a, results1.9b) %>% 
  left_join(results1.9c) %>%
  rbind(vec) 

RQ1.9d <- RQ1.9c
Anova(RQ1.9d, type = 3)
summary(RQ1.9d)
furnival(RQ1.9d)
post5$yhat <- predict(RQ1.9d, level = 0)
bias(post5$yhat, post5$EW.dens)

# LW dens post 1972
RQ1.10a <- lme(
  LW.dens ~ YST * Sample.HT + PDSIavg + LW.dens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post10
)
RQ1.10b <- lme(
  LW.dens ~ YST * Sample.HT + PDSIavg + LW.dens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = post5
)
RQ1.10c <- lme(
  LW.dens ~ YST * Sample.HT + PDSIavg + LW.dens.pre72 + GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  correlation = corCAR1(form = ~ 1),
  data = post5
)
plot(RQ1.10c)

post10$yhat <- predict(RQ1.10a, level = 0)
bias(post10$yhat, post10$LW.dens)

post5$yhat <- predict(RQ1.10b, level = 0)
bias(post5$yhat, post5$LW.dens)

post5$yhat <- predict(RQ1.10c, level = 0)
bias(post5$yhat, post5$LW.dens)

results1.10a <- results.fn(RQ1.10a) %>% 
  rename("OG" = value)
results1.10b <- results.fn(RQ1.10b) %>% 
  rename("OGPost5" = value)
results1.10c <- results.fn(RQ1.10c) %>% 
  rename("Post5CAR" = value)

vec <- c("NA","NA", furnival(RQ1.10a)$coefd[2,1],
         furnival(RQ1.10b)$coefd[2,1],
         furnival(RQ1.10c)$coefd[2,1])

results1.10 <- left_join(results1.10a, results1.10b) %>% 
  left_join(results1.10c) %>%
  rbind(vec)

Anova(RQ1.10d, type = 3)
RQ1.10d <- update(RQ1.10c, LW.dens ~ YST + PDSIavg + LW.dens.pre72)
summary(RQ1.10d)
furnival(RQ1.10d)
post5$yhat <- predict(RQ1.10d, level = 0)
bias(post5$yhat, post5$LW.dens)

# Prediction graphs

# Graph predictions from lme2- problem with predict(). https://stats.stackexchange.com/questions/29513/error-in-getting-predictions-from-a-lme-object
Sample.HT <- unique(df4$Sample.HT)
YST <- seq(0, 9)
others <- df4 %>% 
  select(YST, PDSIavg) %>% 
  slice(1:10)

ringDens10 <- df4 %>% 
  group_by(Sample.HT) %>% 
  summarise(ringDens10 = quantile(ringDens10, probs = .5, na.rm = TRUE))

newdata <- expand.grid(YST, Sample.HT) %>% 
  rename("YST" = Var1, "Sample.HT" = Var2) %>% 
  left_join(others, by = "YST") %>% 
  left_join(ringDens10, by = "Sample.HT") %>% 
  na.omit()

newdata$yhat <- predict(lme1, newdata, level = 0)

# Doesn't pick up on a big spike 4-5 years after thinning
ggplot(newdata) +
  geom_line(aes(x = YST, y = yhat, color = factor(Sample.HT)), 
            linetype = 2) +
  # geom_smooth(aes(x = YST, y = ringDens, color = factor(Sample.HT)), 
  #             se = FALSE, data = filter(df4, YST <= 10)) 
  stat_summary(aes(x = YST, y = ringDens, color = factor(Sample.HT)),
               fun.y = median, geom = "point",
              data = filter(df4, YST <= 9))
  # geom_point(aes(x = YST, y = ringDens, color = factor(Sample.HT)),
  #             se = FALSE, data = filter(df4, YST <= 10))

# New stuff 2/1 -----------------------------------------------------------

df3 <- df2a %>% 
  filter(Sample.HT <= 24)

AddPre10 <- df3 %>% 
  filter(year >= 1953 & year <= 1962) %>% 
  group_by(standtree) %>% 
  summarise(BAI10 = mean(BAI, na.rm = TRUE),
            ringDens10 = mean(ringDens, na.rm = TRUE),
            LW.prop10 = mean(LW.prop, na.rm = TRUE),
            LW.dens10 = mean(LW.dens, na.rm = TRUE),
            EW.dens10 = mean(EW.dens, na.rm = TRUE)) 

df3 <- df3 %>% 
  left_join(AddPre10, by = "standtree")

library(mgcv)
mod1 <- lme(ringDens ~ newring * Sample.HT + PDSIavg,
            random = ~1|TreeID,
            na.action = na.omit,
            data = df3)
furnival(mod1)
Anova(mod1, type = 3)
summary(mod1)

df4 <- df3 %>% 
  filter(year >= 1962) %>% 
  na.omit()

mod2 <- lme(ringDens ~ Sample.HT + prec_3 + GDD_3 + GDD_4 +
              ringDens10,
            random = ~1|TreeID,
            # weights = varPower(form = ~RW),
            # correlation = corCAR1(form = ~1),
            data = df4)
furnival(mod2)
Anova(mod2, type = 3)
summary(mod2)

plot(mod2)
df4$res <- resid(mod2, type = "pearson")
ggplot(df4) +
  geom_smooth(aes(x = RW, y = res))

df4$yhat <- predict(mod2, level = 0)
bias(df4$yhat, df4$ringDens)

#           E        |E|       RMSE         E%       |E|% 
#   1.1504739 36.5222598 47.4803426  0.2022738  8.1680773 

ggplot(df4) +
  # geom_smooth(aes(x = year, y = yhat, color = factor(Sample.HT))) +
  stat_summary(aes(x = year, y = yhat, color = factor(Sample.HT)),
               fun.y = median, geom = "line")

Sample.HT <- unique(df3$Sample.HT)

ringDens10 <- df3 %>% 
  group_by(Sample.HT) %>% 
  summarise(ringDens10 = quantile(ringDens10, probs = .5, na.rm = TRUE))

others <- df4 %>% 
  select(YST, PDSIavg) %>% 
  slice(1:10)



newdata <- expand.grid(YST, Sample.HT) %>% 
  rename("YST" = Var1, "Sample.HT" = Var2) %>% 
  left_join(others, by = "YST") %>% 
  left_join(ringDens10, by = "Sample.HT") %>% 
  na.omit()

newdata$yhat <- predict(lme1, newdata, level = 0)

# Old Research Question 3 -----------------------------------------------------

climWide <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", header = TRUE) 

climWideLag1 <- climWide %>% 
  mutate_at(vars(-grow.yr), lag) %>% 
  rename_at(vars(-grow.yr), funs(paste(.,"lag1", sep = ""))) 

climWideLag2 <- climWideLag1 %>% 
  mutate_at(vars(-grow.yr), lag) %>% 
  rename_at(vars(-grow.yr), funs(paste(.,"lag2", sep = ""))) 

climWideLag <- climWide %>% 
  left_join(climWideLag1, by = "grow.yr") %>% 
  left_join(climWideLag2, by = "grow.yr")

climWide2 <- climWide %>% 
  mutate(PRCPlag1 = lag(PRCPtot, n = 1),
         PRCPlag2 = lag(PRCPtot, n = 2),
         PRCPlag3 = lag(PRCPtot, n = 3),
         PRCPlag4 = lag(PRCPtot, n = 4),
         PRCPlag5 = lag(PRCPtot, n = 5),
         PRCPlag6 = lag(PRCPtot, n = 6),
         PRCPlag7 = lag(PRCPtot, n = 7),
         PRCPlag8 = lag(PRCPtot, n = 8),
         PRCPlag9 = lag(PRCPtot, n = 9),
         PRCPlag10 = lag(PRCPtot, n = 10),
         dummy = rnorm(98, 545, 158))

df3 <- df2a %>% 
  filter(Sample.HT <= 24) %>% 
  select(standtree, SiteID, PlotID, TreeID, Position, Sample.HT, year, newring, 
         BAI, ringDens, LW.prop, EW.dens, LW.dens) %>% 
  left_join(climWide2, by = c("year" = "grow.yr")) %>% 
  mutate(dummy2 = rnorm(15932, 545, 189))

ggplot(filter(df3, PlotID == 3)) +
  geom_line(aes(x = PDSIavg, y = ringDens, color = factor(TreeID))) +
  facet_wrap(~Sample.HT)

testmod <- gam(
  ringDens ~ s(PDSIavg) + Sample.HT + s(newring),
  data = df3
)
summary(testmod)
plot.gam(testmod)

testmod <- lm(
  ringDens ~ dummy,
  data = df3
)
summary(testmod)

testmod <- lme(
  ringDens ~ Sample.HT + newring + GSL,
    random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df3
)
Anova(testmod, type = 3)

df4 <- df3 %>% 
  select(standtree, dummy, dummy2)

# Use df3 for the full dataframe (comment the pre-62s); df4 is subsetted (remember to un-comment the pre62's). Results are really about the same regarding weather.

# Need to define TreeID as a factor for gam to work
str(df2a)
df2a$TreeID <- as.factor(df2a$TreeID)

# df3 <- df2a %>% 
#   filter(Sample.HT <= 24) %>% 
#   select(standtree, SiteID, PlotID, TreeID, Position, Sample.HT, year, newring, 
#          BAI, ringDens, LW.prop, EW.dens, LW.dens) %>% 
#   left_join(climWideLag, by = c("year" = "grow.yr"))

df3 <- df2a %>% 
  filter(Sample.HT <= 24) %>% 
  select(standtree, SiteID, PlotID, TreeID, Position, Sample.HT, year, newring, 
         BAI, ringDens, LW.prop, EW.dens, LW.dens) %>% 
  left_join(climWide2, by = c("year" = "grow.yr")) %>% 
  mutate(dummy2 = rnorm(15932, 545, 189))

# Summary stats below are from full df2a based on purrr with 80/20 method

coef.fn <- function(y,x){
  summary(y)$p.coef[x]
}
pval.fn <- function(y,x){
  summary(y)$p.pv[x]
}

testmod <- gam(
  ringDens ~ Sample.HT + 
    PRCPtot + PRCPlag1 + PRCPlag2 + PRCPlag3 + PRCPlag4 + PRCPlag5 + PRCPlag6 + PRCPlag7 + PRCPlag8 + PRCPlag9 +
    PRCPlag10 + dummy + dummy2 +  
    s(newring, k = 5),
  data = df3
)
summary(testmod)

testmod <- lme(
  ringDens ~ Sample.HT + 
    PRCPtot + PRCPlag1 + PRCPlag2 + PRCPlag3 + PRCPlag4 + PRCPlag5 + PRCPlag6 + PRCPlag7 + PRCPlag8 + PRCPlag9 +
    PRCPlag10 + dummy + dummy2,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df3
)
Anova(testmod, type = 3)

library(mgcv)
RQ3.1 <- gam(
  BAI ~ Sample.HT + prec_1 + prec_2 + prec_3 + prec_4 + 
    temp_1 + temp_2 + temp_3 + temp_4 + s(newring, k = 5) +
    # BAI.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)

RQ3.1 <- gam(
  BAI ~ Sample.HT + 
    PRCP.1lag1lag2 + PRCP.2lag1lag2 + PRCP.3lag1lag2 + PRCP.4lag1lag2 +
    TAVG.1lag1lag2 + TAVG.2lag1lag2 + TAVG.3lag1lag2 + TAVG.4lag1lag2 +
    PRCP.1lag1 + PRCP.2lag1 + PRCP.3lag1 + PRCP.4lag1 +
    TAVG.1lag1 + TAVG.2lag1 + TAVG.3lag1 + TAVG.4lag1 +
    PRCP.1 + PRCP.2 + PRCP.3 + PRCP.4 + 
    TAVG.1 + TAVG.2 + TAVG.3 + TAVG.4 + 
    s(newring, k = 5) +
    s(TreeID, bs = "re"),
  data = df3
)
anova.gam(RQ3.1)
summary(RQ3.1)

RQ3.1 <- gam(
  BAI ~ Sample.HT + scale(prec_1, center = TRUE, scale = TRUE) + 
    scale(prec_2, center = TRUE, scale = TRUE) + 
    scale(prec_3, center = TRUE, scale = TRUE) + 
    scale(prec_4, center = TRUE, scale = TRUE) + 
    scale(temp_1, center = TRUE, scale = TRUE) + 
    scale(temp_2, center = TRUE, scale = TRUE) + 
    scale(temp_3, center = TRUE, scale = TRUE) + 
    scale(temp_4, center = TRUE, scale = TRUE) + 
    s(newring, k = 5) +
    # BAI.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)

summary(RQ3.1)
plot.gam(RQ3.1)

RQ3.1results <- data.frame(covariate = c("prec1", "prec2", "prec3", "prec4", 
                                         "temp1", "temp2", "temp3", "temp4"),
                           rownum = seq(3,10)) %>% 
  mutate(Pval = pval.fn(RQ3.1,rownum),
         Effect = ifelse(Pval < .05, coef.fn(RQ3.1,rownum), 0))

ggplot(RQ3.1results) +
  geom_bar(aes(x = covariate, y = Effect), stat = "identity") +
  ggtitle("BAI response")

RQ3.2 <- gam(
  ringDens ~ Sample.HT + PRCP.1 + PRCP.2 + PRCP.3 + PRCP.4 + TAVG.1 + 
    TAVG.2 + TAVG.3 + TAVG.4 + s(newring, k = 5) +
    # ringDens.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)
summary(RQ3.2)
RQ3.2 <- gam(
  ringDens ~ Sample.HT + scale(PRCP.1, center = TRUE, scale = TRUE) + 
    scale(PRCP.2, center = TRUE, scale = TRUE) + 
    scale(PRCP.3, center = TRUE, scale = TRUE) + 
    scale(PRCP.4, center = TRUE, scale = TRUE) + 
    scale(TAVG.1, center = TRUE, scale = TRUE) + 
    scale(TAVG.2, center = TRUE, scale = TRUE) + 
    scale(TAVG.3, center = TRUE, scale = TRUE) + 
    scale(TAVG.4, center = TRUE, scale = TRUE) + 
    s(newring, k = 5) +
    # BAI.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)
summary(RQ3.2)
plot.gam(RQ3.2)

RQ3.2results <- data.frame(covariate = c("prec1", "prec2", "prec3", "prec4", 
                                         "temp1", "temp2", "temp3", "temp4"),
                           rownum = seq(3,10)) %>% 
  mutate(Pval = pval.fn(RQ3.2,rownum),
         Effect = ifelse(Pval < .05, coef.fn(RQ3.2,rownum), 0))

ggplot(RQ3.2results) +
  geom_bar(aes(x = covariate, y = Effect), stat = "identity") +
  ggtitle("Mean ring density response")

RQ3.3 <- gam(
  LW.prop ~ Sample.HT + prec_1 + prec_2 + prec_3 + prec_4 + 
    temp_1 + temp_2 + temp_3 + temp_4 + s(newring, k = 5) +
    # LW.prop.pre62 + 
    s(TreeID, bs = "re"),
  data = df3
)

RQ3.3 <- gam(
  LW.prop ~ Sample.HT + scale(prec_1, center = TRUE, scale = TRUE) + 
    scale(prec_2, center = TRUE, scale = TRUE) + 
    scale(prec_3, center = TRUE, scale = TRUE) + 
    scale(prec_4, center = TRUE, scale = TRUE) + 
    scale(temp_1, center = TRUE, scale = TRUE) + 
    scale(temp_2, center = TRUE, scale = TRUE) + 
    scale(temp_3, center = TRUE, scale = TRUE) + 
    scale(temp_4, center = TRUE, scale = TRUE) + 
    s(newring, k = 5) +
    # BAI.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)
summary(RQ3.3)
plot.gam(RQ3.3)

RQ3.3results <- data.frame(covariate = c("prec1", "prec2", "prec3", "prec4", 
                                         "temp1", "temp2", "temp3", "temp4"),
                           rownum = seq(3,10)) %>% 
  mutate(Pval = pval.fn(RQ3.3,rownum),
         Effect = ifelse(Pval < .05, coef.fn(RQ3.3,rownum), 0))

ggplot(RQ3.3results) +
  geom_bar(aes(x = covariate, y = Effect), stat = "identity") +
  ggtitle("Latewood proportion response")


RQ3.4 <- gam(
  EW.dens ~ Sample.HT + prec_1 + prec_2 + prec_3 + prec_4 + 
    temp_1 + temp_2 + temp_3 + temp_4 + s(newring, k = 5) +
    # EW.dens.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)

RQ3.4 <- gam(
  EW.dens ~ Sample.HT + scale(prec_1, center = TRUE, scale = TRUE) + 
    scale(prec_2, center = TRUE, scale = TRUE) + 
    scale(prec_3, center = TRUE, scale = TRUE) + 
    scale(prec_4, center = TRUE, scale = TRUE) + 
    scale(temp_1, center = TRUE, scale = TRUE) + 
    scale(temp_2, center = TRUE, scale = TRUE) + 
    scale(temp_3, center = TRUE, scale = TRUE) + 
    scale(temp_4, center = TRUE, scale = TRUE) + 
    s(newring, k = 5) +
    # BAI.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)
summary(RQ3.4)
plot.gam(RQ3.4)

RQ3.4results <- data.frame(covariate = c("prec1", "prec2", "prec3", "prec4", 
                                         "temp1", "temp2", "temp3", "temp4"),
                           rownum = seq(3,10)) %>% 
  mutate(Pval = pval.fn(RQ3.4,rownum),
         Effect = ifelse(Pval < .05, coef.fn(RQ3.4,rownum), 0))

ggplot(RQ3.4results) +
  geom_bar(aes(x = covariate, y = Effect), stat = "identity") +
  ggtitle("Earlywood density response")

RQ3.5 <- gam(
  LW.dens ~ Sample.HT + prec_1 + prec_2 + prec_3 + prec_4 + 
    temp_1 + temp_2 + temp_3 + temp_4 + s(newring, k = 5) +
    # LW.dens.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)

RQ3.5 <- gam(
  LW.dens ~ Sample.HT + scale(prec_1, center = TRUE, scale = TRUE) + 
    scale(prec_2, center = TRUE, scale = TRUE) + 
    scale(prec_3, center = TRUE, scale = TRUE) + 
    scale(prec_4, center = TRUE, scale = TRUE) + 
    scale(temp_1, center = TRUE, scale = TRUE) + 
    scale(temp_2, center = TRUE, scale = TRUE) + 
    scale(temp_3, center = TRUE, scale = TRUE) + 
    scale(temp_4, center = TRUE, scale = TRUE) + 
    s(newring, k = 5) +
    # BAI.pre62 +
    s(TreeID, bs = "re"),
  data = df3
)
summary(RQ3.5)
plot.gam(RQ3.5)

RQ3.5results <- data.frame(covariate = c("prec1", "prec2", "prec3", "prec4", 
                                         "temp1", "temp2", "temp3", "temp4"),
                           rownum = seq(3,10)) %>% 
  mutate(Pval = pval.fn(RQ3.5,rownum),
         Effect = ifelse(Pval < .05, coef.fn(RQ3.5,rownum), 0))

ggplot(RQ3.5results) +
  geom_bar(aes(x = covariate, y = Effect), stat = "identity") +
  ggtitle("Latewood density response")

# Plots -------------------------------------------------------------------

# Average ring width, useful as a double-check to make sure drought years are right
d <- df2 %>% 
  group_by(year) %>% 
  summarise(ringDens = mean(ringDens, na.rm = TRUE))
ggplot(d) +
  geom_bar(aes(x = year, y = ringDens), stat = "identity") +
  geom_vline(xintercept = 2002) +
  geom_vline(xintercept = 1996) +
  geom_vline(xintercept = 1989) +
  geom_vline(xintercept = 1977) +
  geom_vline(xintercept = 1963) +
  geom_vline(xintercept = 1956) +
  geom_vline(xintercept = 1951)

# Ring density and ring width over years at different Positions
ggplot(df2) +
  # geom_point(aes(x=year, y=ringDens), col = "gray50") +
  geom_smooth(aes(x=year, y=ringDens, color = factor(Position)),
              method = "gam",
              formula = y ~ s(x, k = 10),
              # se = FALSE,
              size = 1) 

# Plot of weather variables
ggplot(w.sum) +
  geom_smooth(aes(x = grow.yr, y = rain), method = "loess")

ggplot(precip) +
  # geom_line(aes(x = month, y = snow)) +
  stat_summary(aes(x = grow.yr, y = prec), geom = "bar", fun.y = sum) +
  geom_vline(xintercept = 2002) +
  geom_vline(xintercept = 1996) +
  geom_vline(xintercept = 1989) +
  geom_vline(xintercept = 1977) +
  geom_vline(xintercept = 1974) +
  geom_vline(xintercept = 1971) +
  geom_vline(xintercept = 1963) +
  geom_vline(xintercept = 1956) +
  geom_vline(xintercept = 1951)

# Overall ringDens vs year for all samples
ggplot(df) +
  # geom_point(aes(x=newring, y=ringDens), col = "gray50") +
  geom_smooth(aes(x=year, y=ringDens, color = factor(Position)),
              # method = "gam",
              # se = FALSE,
              size = 1) +
  ylab("Mean ring density (kg/m^3)") +
  xlab("Ring Number") +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size =16))  
+
  facet_wrap(~GSL)

# General trend of min, EW, mean, LW, and max densities
ggplot(df2) +
  geom_smooth(aes(x = newring, y = minDens), 
              method = "gam", formula = y ~ s(x, k = 3)) +
  geom_smooth(aes(x = newring, y = maxDens), 
              method = "gam", formula = y ~ s(x, k = 3)) +
  geom_smooth(aes(x = newring, y = EW.dens), 
              method = "gam", formula = y ~ s(x, k = 3)) +
  geom_smooth(aes(x = newring, y = LW.dens), 
              method = "gam", formula = y ~ s(x, k = 3)) +
  geom_smooth(aes(x = newring, y = ringDens), 
              method = "gam", formula = y ~ s(x, k = 3)) +
  ylab("Ring Density (min, EW.avg, average, LW.avg, max") +
  xlab("newring") 

# Look at effect of crown distance. Change filter to look at effect of year, newring, and mm range
ggplot(filter(df2, mmRingStart >= 50 & mmRingStart < 60)) +
  geom_smooth(aes(x = RelHT, y = ringDens),
              method = "gam", formula = y ~ s(x, k = 5))


# Tree ring with height profile -------------------------------------------

dfmmRing <- df2 %>% 
  filter(TreeID == "1-146")

# 12-21 one of the best
dfmmRing <- sumProfile2 %>% 
  filter(TreeID == "12-21")

# ggplot(dfmmRing) +
#   geom_point(aes(x = factor(Sample.HT), y = mmRingStart, color = year)) +
#   coord_flip()

ggplot(dfmmRing) +
  geom_line(aes(x = Sample.HT, 
                y = mmRingStart, 
                color = year, 
                group = year)) +
  coord_flip() +
  geom_vline(xintercept = dfmmRing$CBH)

# Add label to every 10th year- %% means remainder
dfmmRing2 <- dfmmRing %>% 
  group_by(year) %>% 
  filter(year%%10 == 0 & mmRingStart == max(mmRingStart))
  
ggplot(dfmmRing, aes(x = Sample.HT, y = mmRingStart)) +
  geom_line(aes(color = year, 
                group = year)) +
  coord_flip() +
  geom_label(data = dfmmRing2, aes(label = year))

# Labels not that helpful. Better to change color of certain years- 1962, other thinnings. Flip color palette

dfmmRing <- dfmmRing %>% 
  mutate(thinyr = ifelse(year == 1962 | year == 1972 | 
                           year == 1982 | year == 1992 | 
                           year ==2002, "Y", "N"))

ggplot(dfmmRing) +
  geom_line(aes(x = Sample.HT, 
                y = mmRingStart, 
                color = thinyr,
                group = year)) +
  coord_flip() +
  geom_vline(xintercept = dfmmRing$CBH)


# Multi-plots -------------------------------------------------

df2 <- df2a %>% 
  filter(Sample.HT <= 24)

climWide <- read.csv("AnalysisReadyData/climWide.csv", header = TRUE)

# The one Dave was looking for- with climate at top
precJuly <- read.csv("AnalysisReadyData/precJuly.csv", header = TRUE)
climWide <- left_join(climWide, precJuly, by = "grow.yr")

ggplot(climWide) +
  geom_line(aes(x = grow.yr, y = PDSIavg))
ggplot(climWide) +
  geom_line(aes(x = grow.yr, y = precJuly))
ggplot(climWide) +
  geom_line(aes(x = grow.yr, y = prectot))

p1 <- ggplot(filter(climWide, grow.yr >= 1940 & grow.yr <= 2016)) +
  geom_line(aes(x = grow.yr, y = PDSIavg)) +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
                            linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("PDSI") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p2 <- ggplot(filter(df2, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = ringDens), fun.y = mean, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("Mean ring density") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p3 <- ggplot(filter(df2, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = LW.prop), fun.y = mean, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("Latewood proportion") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p4 <- ggplot(filter(df2, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = EW.dens), fun.y = mean, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("Mean EW density") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p5 <- ggplot(filter(df2, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = LW.dens), fun.y = mean, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  ylab("Mean LW density") +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)
g5 <- ggplotGrob(p5)
g <- rbind(g1, g2, g3, g4, g5, size = "first")
grid.newpage()
grid.draw(g)


# Similar but with newring

p1 <- ggplot(filter(df2, newring >= 10 & Sample.HT <= 24)) +
  stat_summary(aes(x = newring, y = ringDens, color = factor(Sample.HT)),
               fun.y = mean, geom = "line") +
  # geom_smooth(aes(x = newring, y = ringDens, color = factor(Sample.HT)),
  #              formula = y ~ s(x, k = 3)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p2 <- ggplot(filter(df2, newring >= 10 & Sample.HT <= 24)) +
  stat_summary(aes(x = newring, y = LW.prop, color = factor(Sample.HT)),
               fun.y = mean, geom = "line") +
  # geom_smooth(aes(x = newring, y = LW.prop, color = factor(Sample.HT)),
  #             formula = y ~ s(x, k = 3)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p3 <- ggplot(filter(df2, newring >= 10 & Sample.HT <= 24)) +
  stat_summary(aes(x = newring, y = EW.dens, color = factor(Sample.HT)),
               fun.y = mean, geom = "line") +
  # geom_smooth(aes(x = newring, y = EW.dens, color = factor(Sample.HT)),
  #             formula = y ~ s(x, k = 3)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p4 <- ggplot(filter(df2, newring >= 10 & Sample.HT <= 24)) +
  stat_summary(aes(x = newring, y = LW.dens, color = factor(Sample.HT)),
               fun.y = mean, geom = "line") +
  # geom_smooth(aes(x = newring, y = LW.dens, color = factor(Sample.HT)),
  #             formula = y ~ s(x, k = 3)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)
g <- rbind(g1, g2, g3, g4, size = "first")
grid.newpage()
grid.draw(g)

# Models ------------------------------------------------------------------

library(lme4)
library(lmerTest)
library(mgcv)

YST <- data.frame(year = seq(1963, 2016), 
                  YST = c(rep(seq(1,10), 4), seq(1,14)))

library(nlme)
library(mgcv)

df.62 <- df2 %>% 
  filter(year >= 1962 & year <= 2014) %>% 
  filter(newring >= 20) %>% 
  left_join(YST, by = "year") %>% 
  na.omit()

lm.list <- nlme::lmList(ringDens ~ RelHT + PDSIavg + YST | TreeID, 
                        data = df.62, na.action = na.omit)
lm.list2 <- nlme::lmList(ringDens ~ RelHT + PDSIavg | factor(YST), 
                        data = df.62, na.action = na.omit)

# interval plots are very useful
coef(lm.list)
plot(intervals(lm.list2))

# The model below is good- high r2. But what are the implications of having sample.HT within treeID as random effect? If you just do TreeID, r2 is much lower.
lme1 <- lme(ringDens ~ PDSIavg + YST * factor(GSL), 
            random = ~1|Sample.HT/TreeID,
            data = df.62)
anova(lme1)
df.62$yhat <- predict(lme1, re.form = NA)
modtest <- lm(yhat ~ ringDens, data = df.62)
summary(modtest)

lm1 <- lm(LW.prop ~ YST * factor(GSL) * Sample.HT, data = df.62)
summary(lm1)
anova(lm1)

# Newring basically represents the low frequency variation over the disk's life. Including a newring smoother would remove the thinning and TRP effect. Analyzing the smoothed trend would tell you about the thinning and TRP effect. So, use a GAM instead of this model.
lm2 <- lm(LW.prop ~ newring * Sample.HT + PDSIavg, data = df.62)
summary(lm2)

gam1 <- mgcv::gam(LW.prop ~ s(newring, k = 5) + Sample.HT + PDSIavg, 
                  data = df2)
summary(gam1)
plot.gam(gam1)
# can do other plots than resid vs. fitted...
plot(lme1, form = resid(., type = "p") ~ year)

lme2 <- update(lme1, random = ~1|PlotID/TreeID)
summary(lme2)
intervals(lme2)

anova(lme1, lme2)
# PlotID variance has very wide intervals. Also the anova shows that inclusion does not include the fit. So I'm dropping PlotID random and trying as a fixed effect

lme3 <- lme(ringDens ~year + RW + factor(PlotID) + factor(Sample.HT), 
            random = ~1|TreeID,
            data = df2)
summary(lme3)
anova(lme3)
emmeans(lme3, ~PlotID)

# SHows that PlotID is just barely not significant. Only Plot 10 is significantly different from any others- not enough to justify its inclusion

lme4 <- lme(ringDens ~year + RW + factor(GSL) + factor(Sample.HT), 
            random = ~1|TreeID,
            data = df2)
summary(lme4)
anova(lme4)
emmeans(lme4, ~GSL*Sample.HT)

# Interaction is significant... possible to incorporate into nonlinear models?
lme3 <- lme(ringDens ~ RW + newring * factor(Sample.HT),
            random = ~1|TreeID,
            data = df2)
summary(lme3)
anova(lme3)




# Linear mixed with GSL. GSL is not significant
lme1 <- lmer(ringDens ~ Sample.HT + RW + YST + GSL + PDSIavg + (1|TreeID),
             data = df.62,
             na.action = na.omit)
summary(lme1)

df.62$yhat <- predict(lme1, re.form = NA)

modtest <- lm(yhat ~ ringDens, data = df.62)
summary(modtest)

# Try Plot BA. Plot BA, as a continuous variable, is not significant
PlotBA <- read.csv("IntermediateData/TW_PlotBA.csv", header = TRUE) %>% 
  select(-X)
PlotBA$PlotID <- as.factor(PlotBA$PlotID)

df.62 <- left_join(df.62, PlotBA, by = "PlotID")

lme2 <- lmer(ringDens ~ HT + RW + year + pre62 + BA + (1|TreeID),
             data = df.62,
             na.action = na.omit)
summary(lme2)

# Try PlotID instead of BA or GSL. This is significant, showing that plot has an effect that is not related to density. Likely microsite conditions
str(df2)
lme3 <- lmer(ringDens ~ HT + rain + snow + newring + (1|TreeID),
             data = df,
             na.action = na.omit)
anova(lme3)
summary(lme3)
coef(lme3)

df$yhat <- predict(lme3, re.form = NA)

modtest <- lm(yhat ~ ringDens, data = df)
summary(modtest)

# Add crown ratio. Not significant, but a linear model without the Tree random effect has it as significant. Conclusion: Can't include terms in the model that occur at the same level as the random effect. This true? If so, it it ok to have plot level variables?

TW.tree <- read.csv("IntermediateData/TreeLevelDataTW.csv", header = TRUE) %>% 
  mutate(TreeID = paste(PlotID, TreeID, sep = "-"),
         CR = (HT - CBH)/HT) %>% 
  select(TreeID, CR)

df.62.cr <- left_join(df.62, TW.tree, by = "TreeID")

lme4 <- lmer(ringDens ~ HT + RW + year + pre62 + PlotID + CR + (1|TreeID),
             data = df.62.cr,
             na.action = na.omit)
anova(lme4)
summary(lme4)

lm3 <- lm(ringDens ~ HT + RW + year + pre62 + PlotID + CR,
            data = df.62.cr,
            na.action = na.omit)
anova(lm3)

# Abandon the tree level stuff. Check letting others vary by TreeID. Letting RW vary with TreeID makes year insignificant

lme5 <- lmer(ringDens ~ HT + RW + year + pre62 + PlotID + 
               (1|TreeID) + (RW|TreeID),
             data = df.62,
             na.action = na.omit)
anova(lme5)
summary(lme5)
coef(lme5)

# Best model?
lme6 <- lmer(ringDens ~ RelHT + year + PDSIavg + (1|TreeID),
             data = df2,
             na.action = na.omit)
summary(lme6)
anova(lme6)

df2$yhat <- predict(lme6, newdata = df2, re.form = NA)
RMSE <- sqrt(mean((df2$yhat - df2$ringDens)^2, na.rm = TRUE))

modtest <- lm(yhat ~ ringDens, data = df.62)
summary(modtest)

newyear <- seq(from = 1918, to = 2018, by = 1)
# newRelHT <- quantile(df2$RelHT, probs = c(.25, .5, .75))
newRelHT <- c(.5, 1, 1.5)

newdat <- expand.grid(newyear, newRelHT) %>%
  rename(year = Var1, RelHT = Var2)

newdat2 <- df2 %>% 
  group_by(year) %>%
  summarise(PDSIavg = unique(PDSIavg)) %>% 
  right_join(newdat)

newdat2$yhat <- predict(lme6, newdata = newdat2, re.form = NA)

ggplot(newdat2) +
  geom_line(aes(x = year, y = yhat, color = factor(RelHT)), size = 1) 
  # geom_point(aes(x = year, y = yearDens), data = twBendFull, color = 'gray50') +
  # theme(legend.position = "none") 

  # facet_wrap(~GSLmetric, nrow = 2) +
  # xlab("year Number") +
  # ylab("Predicted yearDens (GPa)") +
  # theme(axis.title.x = element_text(size = 18)) +
  # theme(axis.title.y = element_text(size = 18)) +
  # theme(axis.text.x = element_text(size = 14)) +
  # theme(axis.text.y = element_text(size = 14))

m1 <- gam(ringDens ~ year + s(RelHT), data = df2)
plot.gam(m1)

# Density vs relht and year

ggplot(filter(df2, TreeID == "14-212" & year == 1985)) + 
  geom_smooth(aes(x = Sample.HT, y = RW))
  


library(mgcv)

gam1 <- gam(ringDens ~ s(year, k = 5), data = df2)
summary(gam1)
plot.gam(gam1)

gam2 <- gam(ringDens ~ s(RelHT, k = 5) + s(year, k = 5), data = df2)
summary(gam2)

gam3 <- gam(ringDens ~ s(RelHT, k = 5) + s(year, k = 5) + s(factor(df2$TreeID), bs = "re"), data = df2)
summary(gam3)

gam4 <- gam(ringDens ~ s(RelHT, k = 5) + s(year, k = 5) + 
              s(factor(df2$TreeID), bs = "re") + s(BAI, k = 5), data = df2)
summary(gam4)

gam5 <- gam(ringDens ~ RelHT + s(year, k = 10) + s(factor(df2$TreeID), bs = "re") + BAI,
            data = df2)
summary(gam5)
gam.check(gam5)
anova(gam5)

# Predicted value over year- not residual
df2$yhat <- predict(gam5, exclude = "s(TreeID)")
ggplot(df2) +
  geom_smooth(aes(x = year, y = yhat)) 

# Detrend -----------------------------------------------------------------

# 1: detrend with dplR to produce index
# 2: Remove climate variability

library(dplR)
df2 <- df %>% 
  filter(Position == 45)

ggplot(df2) +
  # geom_point(aes(x = year, y = RW, color = factor(GSL))) +
  stat_summary(aes(x = year, y = RW),  
               fun.y = mean, geom = "line") +
  theme(legend.position = "none")

df3 <- df2 %>%
  ungroup() %>% 
  data.frame() %>% 
  select(TreeID,year,Position,RW) %>%
  unite(NewID, TreeID, Position, remove = TRUE) %>%
  spread(key = NewID, value = RW) 

rownames(df3) <- df3$year

df3.detrend <- detrend(df3, method = "Spline")
df3.detrend <- detrend(df3, method = "ModNegExp", verbose = TRUE)

df3.detrend.long <- df3.detrend %>% 
  mutate(year = as.numeric(row.names(df3.detrend))) %>% 
  select(year, everything()) %>% 
  gather("standtree", "RWI", 2:51) %>% 
  na.omit() %>% 
  group_by(standtree) %>% 
  mutate(cambialAge = year - min(year) + 1)

ggplot(df3.detrend.long) +
  stat_summary(aes(x = year, y = RWI),  
               fun.y = mean, geom = "line") +
  theme(legend.position="none")  

climWide <- read.csv("AnalysisReadyData/climWide.csv", header = TRUE)
df4 <- left_join(df3.detrend.long, climWide, by = c("year" = "grow.yr"))

df5 <- df4 %>% 
  mutate(newRWI = RWI - (.0007214 * prectot) - (-0.0990936 * temp_4),
         newRWI2 = newRWI - (mean(newRWI, na.rm = TRUE) - 1))

mean(df5$newRWI2, na.rm = TRUE)

ggplot(df5) +
  geom_point(aes(x = year, y = newRWI2, color = factor(standtree))) +
  theme(legend.position="none")

ggplot(df5) +
  stat_summary(aes(x = year, y = newRWI2),  
               fun.y = mean, geom = "line", color = "red") +
  stat_summary(aes(x = year, y = RWI),  
               fun.y = mean, geom = "line", color = "blue")

ggplot(df5) +
  geom_smooth(aes(x = year, y = RWI), color = "blue", se = FALSE) +
  geom_smooth(aes(x = year, y = newRWI2), color = "red", se = FALSE)
# Messin around -----------------------------------------------------------


library(MASS)
boxcox(BAI ~ Sample.HT * YST + PDSIavg + BAI.pre62 + GSL, data = df4, lambda = seq(-.25,.25, length = 10))
library(boxcoxmix)
optim.boxcox(BAI ~ Sample.HT * YST + PDSIavg + BAI.pre62 + GSL, 
             data = df4,
             groups = 53, 
             steps = 500, 
             tol = 0.5,
             start = "gq", 
             EMdev.change = 1e-04, 
             find.in.range = c(-3, 3), 
             s = 60,
             plot.opt = 3, 
             verbose = TRUE, 
             noformat = FALSE)

ARcheck <- df3 %>%
  group_by(TreeID, Sample.HT) %>% 
  select(standtree, TreeID, Sample.HT, BAI, ringDens, 
         LW.prop, EW.dens, LW.dens) %>% 
  mutate(lagBAI = lag(BAI, 1),
         lagringDens = lag(ringDens, 1),
         lagLWprop = lag(LW.prop, 1),
         lagEWdens = lag(EW.dens, 1),
         lagLWdens = lag(LW.dens, 1)) %>% 
  na.omit() %>% 
  nest(-standtree)

cor.BAI <- function(x){
  cor(x$BAI, x$lagBAI)
}
cor.ringDens <- function(x){
  cor(x$ringDens, x$lagringDens)
}
cor.LW.prop <- function(x){
  cor(x$LW.prop, x$lagLWprop)
}
cor.EW.dens <- function(x){
  cor(x$EW.dens, x$lagEWdens)
}
cor.LW.dens <- function(x){
  cor(x$LW.dens, x$lagLWdens)
}

ARcheck2 <- ARcheck %>% 
  mutate(BAIcor = map(data, cor.BAI),
         ringDenscor = map(data, cor.ringDens),
         LWpropcor = map(data, cor.LW.prop),
         EWdenscor = map(data, cor.EW.dens),
         LWdenscor = map(data, cor.LW.dens)) %>% 
  unnest(BAIcor, ringDenscor, LWpropcor, EWdenscor, LWdenscor)
hist(ARcheck2$BAIcor)
hist(ARcheck2$ringDenscor)
hist(ARcheck2$LWpropcor)
hist(ARcheck2$EWdenscor)
hist(ARcheck2$LWdenscor)
getVarCov(RQ1.2)

%>% 
  mutate(x = BAI - mean(BAI),
         y = lag1 - mean(lag1),
         xy = x * y,
         x2 = x^2,
         y2 = y^2)


cor(ARcheck$BAI, ARcheck$lag1)

cor(df3$ringDens, df3$LW.dens)


# Archived stuff from treeclim- use if need to separate by GSL
results30 <- rbind(treeclim.fn(30,.5), 
                   treeclim.fn(30,4.5),
                   treeclim.fn(30,8),
                   treeclim.fn(30,16),
                   treeclim.fn(30,24))

results60 <- rbind(treeclim.fn(60,.5), 
                   treeclim.fn(60,4.5),
                   treeclim.fn(60,8),
                   treeclim.fn(60,16),
                   treeclim.fn(60,24))

results80 <- rbind(treeclim.fn(80,.5), 
                   treeclim.fn(80,4.5),
                   treeclim.fn(80,8),
                   treeclim.fn(80,16),
                   treeclim.fn(80,24))

results100 <- rbind(treeclim.fn(100,.5), 
                    treeclim.fn(100,4.5),
                    treeclim.fn(100,8),
                    treeclim.fn(100,16),
                    treeclim.fn(100,24))

results120 <- rbind(treeclim.fn(120,.5), 
                    treeclim.fn(120,4.5),
                    treeclim.fn(120,8),
                    treeclim.fn(120,16),
                    treeclim.fn(120,24))

results150 <- rbind(treeclim.fn(150,.5), 
                    treeclim.fn(150,4.5),
                    treeclim.fn(150,8),
                    treeclim.fn(150,16),
                    treeclim.fn(150,24))

ggplot(results30, aes(x = month, y = coef, ymin = ci_lower, 
                      ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  ggtitle("BA 30")

ggplot(results60, aes(x = month, y = coef, ymin = ci_lower, 
                      ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  ggtitle("BA 60")

ggplot(results80, aes(x = month, y = coef, ymin = ci_lower, 
                      ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  ggtitle("BA 80")

ggplot(results100, aes(x = month, y = coef, ymin = ci_lower, 
                       ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  ggtitle("BA 100")

ggplot(results120, aes(x = month, y = coef, ymin = ci_lower, 
                       ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  ggtitle("BA 120")

ggplot(results150, aes(x = month, y = coef, ymin = ci_lower, 
                       ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  ggtitle("BA 150")


# Separate functions for diff responses at multiple tree heights

treeclim.BAI <- function(y){
  
  prechron.BAI <- df2 %>%
    filter(Sample.HT == y) %>% 
    select(standtree, year, BAI) %>%
    spread(key = standtree, value = BAI) %>%
    column_to_rownames(var = "year")
  
  prechron.BAI.rwi <- detrend(prechron.BAI, method = "Spline")
  
  chron.BAI <- chron(prechron.BAI.rwi, prewhiten = TRUE) %>%
    rename(BAIIndex = xxxstd)
  
  climtest <- dcc(
    chrono = chron.BAI,
    climate = list(PRCP, TAVG),
    var_names = c("PRCP", "TAVG"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9),
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = y) %>% 
    mutate(quarter = c("ond", "JFM", "AMJ", "JAS", "Annual"))
  
}

resultsBAI <- rbind(
  treeclim.BAI(.5),
  treeclim.BAI(4.5),
  treeclim.BAI(8),
  treeclim.BAI(16),
  treeclim.BAI(24)
) 

BAIplot <- ggplot(resultsBAI, aes(x = quarter, y = coef, ymin = ci_lower, 
                                  ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = resultsBAI$quarter[1:5]) +
  ggtitle("BAI") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Height")

# ringDens

treeclim.ringDens <- function(y){
  
  missYearAvg <- df2 %>%
    filter(Sample.HT == y) %>%
    group_by(year) %>%
    summarise(avgDens = mean(ringDens, na.rm = T)) %>%
    filter(year == 1951 | year == 1956 | year == 1996 | year == 2002)
  
  prechron.ringDens <- df2 %>%
    filter(Sample.HT == y) %>%
    select(standtree, year, ringDens) %>%
    replace_na(list(ringDens = mean(missYearAvg$avgDens))) %>%
    spread(key = standtree, value = ringDens) %>%
    column_to_rownames(var = "year")
  
  prechron.ringDens.rwi <-
    detrend(prechron.ringDens, method = "Spline")
  
  chron.ringDens <- chron(prechron.ringDens.rwi, prewhiten = TRUE) %>%
    rename(DensIndex = xxxstd)
  
  climtest <- dcc(
    chrono = chron.ringDens,
    climate = list(PRCP),
    var_names = c("PRCP"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9),
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = y) %>% 
    mutate(quarter = c("ond", "JFM", "AMJ", "JAS", "Annual"))
}

ringDensResults <- rbind(
  treeclim.ringDens(.5),
  treeclim.ringDens(4.5),
  treeclim.ringDens(8),
  treeclim.ringDens(16),
  treeclim.ringDens(24)
)

ringDensplot <- ggplot(ringDensResults, aes(x = quarter, y = coef, ymin = ci_lower, 
                                            ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = ringDensResults$quarter[1:5]) +
  ggtitle("Ring Density") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Height")
ringDensplot

# LW prop

treeclim.LWprop <- function(y){
  
  missYearAvg <- df2 %>%
    filter(Sample.HT == y) %>%
    group_by(year) %>%
    summarise(avgLWprop = mean(LW.prop, na.rm = T)) %>%
    filter(year == 1951 | year == 1956 | year == 1996 | year == 2002)
  
  prechron.LWprop <- df2 %>%
    filter(Sample.HT == y) %>% 
    select(standtree, year, LW.prop) %>%
    replace_na(list(LW.prop = mean(missYearAvg$avgLWprop))) %>%
    spread(key = standtree, value = LW.prop) %>%
    column_to_rownames(var = "year")
  
  prechron.LWprop.rwi <- detrend(prechron.LWprop, method = "Spline")
  
  chron.LWprop <- chron(prechron.LWprop.rwi, prewhiten = TRUE) %>%
    rename(LWpropIndex = xxxstd)
  
  climtest <- dcc(
    chrono = chron.LWprop,
    climate = list(PRCP, TAVG),
    var_names = c("PRCP", "TAVG"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9),
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = y) %>% 
    mutate(quarter = c("ond", "JFM", "AMJ", "JAS", "Annual"))
}

LWpropResults <- rbind(
  treeclim.LWprop(.5),
  treeclim.LWprop(4.5),
  treeclim.LWprop(8),
  treeclim.LWprop(16),
  treeclim.LWprop(24)
)

LWpropplot <- ggplot(LWpropResults, aes(x = quarter, y = coef, ymin = ci_lower, 
                                        ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = LWpropResults$quarter[1:5]) +
  ggtitle("Latewood proportion") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Height")
LWpropplot

# EW dens

treeclim.EWdens <- function(y){
  
  missYearAvg <- df2 %>%
    filter(Sample.HT == y) %>%
    group_by(year) %>%
    summarise(EWdens = mean(EW.dens, na.rm = T)) %>%
    filter(year == 1951 | year == 1956 | year == 1996 | year == 2002)
  
  prechron.EWdens <- df2 %>%
    filter(Sample.HT == y) %>% 
    select(standtree, year, EW.dens) %>%
    replace_na(list(EW.dens = mean(missYearAvg$EWdens))) %>%
    spread(key = standtree, value = EW.dens) %>%
    column_to_rownames(var = "year")
  
  prechron.EWdens.rwi <- detrend(prechron.EWdens, method = "Spline")
  
  chron.EWdens <- chron(prechron.EWdens.rwi, prewhiten = T) %>%
    rename(EWdensIndex = xxxstd)
  
  climtest <- dcc(
    chrono = chron.EWdens,
    climate = list(PRCP, TAVG),
    var_names = c("PRCP", "TAVG"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9),
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = y) %>% 
    mutate(quarter = c("ond", "JFM", "AMJ", "JAS", "Annual"))
}

EWdensResults <- rbind(
  treeclim.EWdens(.5),
  treeclim.EWdens(4.5),
  treeclim.EWdens(8),
  treeclim.EWdens(16),
  treeclim.EWdens(24)
)

EWdensplot <- ggplot(EWdensResults, aes(x = quarter, y = coef, ymin = ci_lower, 
                                        ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = EWdensResults$quarter[1:5]) +
  ggtitle("Earlywood density") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Height")
EWdensplot

# LW dens

treeclim.LWdens <- function(y){
  
  missYearAvg <- df2 %>%
    filter(Sample.HT == y) %>%
    group_by(year) %>%
    summarise(LWdens = mean(LW.dens, na.rm = T)) %>%
    filter(year == 1951 | year == 1956 | year == 1996 | year == 2002)
  
  prechron.LWdens <- df2 %>%
    filter(Sample.HT == y) %>% 
    select(standtree, year, LW.dens) %>%
    replace_na(list(LW.dens = mean(missYearAvg$LWdens))) %>%
    spread(key = standtree, value = LW.dens) %>%
    column_to_rownames(var = "year")
  
  prechron.LWdens.rwi <- detrend(prechron.LWdens, method = "Spline")
  
  chron.LWdens <- chron(prechron.LWdens.rwi, prewhiten = T) %>%
    rename(LWdensIndex = xxxstd)
  
  climtest <- dcc(
    chrono = chron.LWdens,
    climate = list(PRCP, TAVG),
    var_names = c("PRCP", "TAVG"),
    selection = .sum("PRCP",-10:-12) +
      .sum("PRCP", 1:3) +
      .sum("PRCP", 4:6) +
      .sum("PRCP", 7:9) +
      .sum("PRCP", -10:9),
    boot = "std"
  )
  
  results <- climtest$coef %>% 
    remove_rownames %>% 
    select(-id) %>% 
    mutate(Sample.HT = y) %>% 
    mutate(quarter = c("ond", "JFM", "AMJ", "JAS", "Annual"))
}

LWdensResults <- rbind(
  treeclim.LWdens(.5),
  treeclim.LWdens(4.5),
  treeclim.LWdens(8),
  treeclim.LWdens(16),
  treeclim.LWdens(24)
)

LWdensplot <- ggplot(LWdensResults, aes(x = quarter, y = coef, ymin = ci_lower, 
                                        ymax = ci_upper, color = factor(Sample.HT))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = LWdensResults$quarter[1:5]) +
  ggtitle("Latewood density") +
  ylab("Correlation coefficient") +
  theme(axis.title.x = element_blank()) +
  labs(color = "Height")
LWdensplot
# Archive -----------------------------------------------------------------

# Old graphical summaries
ggplot(filter(df2, GSL == 60 & newring > 2)) +
  stat_summary(aes(x = year, y = LW.prop, color = factor(TreeID)),  
               fun.y = mean, geom = "line") +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003)) 
+
  ggtitle("BAI vs year by GSL, averaged over sample height")

ggplot(df2) +
  stat_summary(aes(x = year, y = BAI, color = factor(Sample.HT)),  
               fun.y = mean, geom = "line") +
  facet_wrap(~GSL) +
  ggtitle("BAI vs year by sample height, GSL facet")

ggplot(df2) +
  stat_summary(aes(x = year, y = BAI, color = factor(GSL)),  
               fun.y = mean, geom = "line") +
  facet_wrap(~Sample.HT) +
  theme_bw() 

ggplot(filter(df2, year >= 1945)) +
  stat_summary(aes(x = year, y = LW.prop, color = factor(GSL)), 
               fun.y = mean, geom = "line") +
  ggtitle("latewood proportion vs year by GSL, averaged over sample position")

ggplot(filter(df2, year >= 1945)) +
  stat_summary(aes(x = year, y = LW.prop, color = factor(Sample.pos)),  
               fun.y = mean, geom = "line") +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003)) +
  ggtitle("Latewood proportion over year, by sample height")

# RQ1 - moving average approach instead of GAM

library(TTY)

BAI.fn <- function(x){SMA(x$BAI, n = 10)}
ringDens.fn <- function(x){SMA(x$ringDens, n = 10)}
LW.prop.fn <- function(x){SMA(x$LW.prop, n = 10)}
EW.dens.fn <- function(x){SMA(x$EW.dens, n = 10)}
LW.dens.fn <- function(x){SMA(x$LW.dens, n = 10)}

df4branch <- df4 %>%
  unnest(data,BAImodel,ringDensmodel,LW.propmodel,EW.densmodel,LW.densmodel) %>% 
  na.omit()

df4branch2 <- df4branch %>%
  filter(year > 1962)

RQ2.1b <- lme(
  BAImodel ~ year + Sample.HT * GSL,
  random = ~ 1 | TreeID,
  data = df4branch2
)
Anova(RQ2.1b, type = 3)
furnival(RQ2.1b)
df4branch2$yhat <- predict(RQ2.1b, level = 0)
bias(df4branch2$yhat, df4branch2$BAI)

RQ2.2b <- lme(
  ringDensmodel ~ year + Sample.HT,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df4branch2
)
Anova(RQ2.2b, type = 3)
furnival(RQ2.2b)
df4branch2$yhat <- predict(RQ2.2b, level = 0)
bias(df4branch2$yhat, df4branch2$ringDensmodel)

RQ2.3b <- lme(
  LW.propmodel ~ year + Sample.HT * GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df4branch2
)
Anova(RQ2.3b, type = 3)
furnival(RQ2.3b)
df4branch2$yhat <- predict(RQ2.3b, level = 0)
bias(df4branch2$yhat, df4branch2$LW.propmodel)

RQ2.4b <- lme(
  EW.densmodel ~ year + Sample.HT * GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df4branch2
)
Anova(RQ2.4b, type = 3)
furnival(RQ2.4b)
df4branch2$yhat <- predict(RQ2.4b, level = 0)
bias(df4branch2$yhat, df4branch2$EW.densmodel)

RQ2.5b <- lme(
  LW.densmodel ~ year + Sample.HT * GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df4branch2
)
Anova(RQ2.5b, type = 3)
furnival(RQ2.5b)
df4branch2$yhat <- predict(RQ2.5b, level = 0)
bias(df4branch2$yhat, df4branch2$LW.densmodel)

# Old Fig 4- good plot!
Fig4a <- ggplot(dfsum, aes(x = factor(GSL), y = BAImean)) +
  geom_boxplot(alpha = 0, size = .75, notch = F, width = .6) +
  geom_jitter(position = position_jitter(0.25), size = 2, alpha = .5, 
              aes(color = factor(Sample.HT))) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 20, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 16)) +
  theme(legend.position = "none") +
  ylab(expression(paste("BAI ", "(cm" ^ "2", ")")))

Fig4b <- ggplot(dfsum, aes(x = factor(GSL), y = RDmean)) +
  geom_boxplot(alpha = 0, size = .75, notch = F, width = .6) +
  geom_jitter(position = position_jitter(0.25), size = 2, alpha = .5, 
              aes(color = factor(Sample.HT))) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 20, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 16)) +
  theme(legend.position = "none") +
  ylab(expression(paste("RD ", "(kg m" ^ "-3", ")")))

Fig4c <- ggplot(dfsum, aes(x = factor(GSL), y = LWPmean)) +
  geom_boxplot(alpha = 0, size = .75, notch = F, width = .6) +
  geom_jitter(position = position_jitter(0.25), size = 2, alpha = .5, 
              aes(color = factor(Sample.HT))) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(size = 20, margin = margin(r = 15)),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 16)) +
  theme(legend.position = "none") +
  xlab(expression(paste("GSL ", "(m" ^ "2 ","ha" ^ "-1", ")"))) +
  ylab("LWP")

Fig4aGrob <- ggplotGrob(Fig4a)
Fig4bGrob <- ggplotGrob(Fig4b)
Fig4cGrob <- ggplotGrob(Fig4c)
Fig4Grob <- rbind(Fig4aGrob, Fig4bGrob, Fig4cGrob, size = "first")
Fig4 <- grid.arrange(Fig4Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/GSL_comparison.tiff", plot = Fig4,
       device = "tiff", 
       width = 32, height = 20, units = "cm")

# Statistics- read Wigley 1984 paper. Some are presented in Pompa-Garcia

chronstats.fn <- function(x,y){
  
  prechron <- df5 %>%
    filter(Sample.HT == x) %>% 
    select(standtree, year, response = y) %>%
    spread(key = standtree, value = response) %>%
    column_to_rownames(var = "year")
  
  prechron.rwi <- detrend(prechron, method = "Spline")
  
  mean.ar1 <- mean(rwl.stats(prechron)$ar1)
  sd.ar1 <- sd(rwl.stats(prechron)$ar1)
  
  results <- rwi.stats(prechron.rwi, prewhiten = T, window.length = 20,
                       window.overlap = 10) %>% 
    select(n.cores, n.tot, rbar.tot, eps, snr) %>% 
    cbind(mean.ar1) %>% 
    cbind(sd.ar1) 
}

chronstats <- rbind(
  chronstats.fn(4.5,32) %>% mutate(response = "BAI"),
  chronstats.fn(4.5,9) %>% mutate(response = "RD"),
  chronstats.fn(4.5,31) %>% mutate(response = "LWP"),
  chronstats.fn(4.5,12) %>% mutate(response = "EWD"),
  chronstats.fn(4.5,11) %>% mutate(response = "MXD")
) %>% 
  select(response, everything())

write.csv(chronstats, "AnalysisReadyData/chronstats.csv", row.names = F)

# GAM method for RQ1

df2$TreeID <- as.factor(df2$TreeID)

df3 <-
  df2 %>%
  select(standtree, PlotID, TreeID, Sample.HT, newring, year, 
         BAI, ringDens, LW.prop, EW.dens, LW.dens, maxDens, GSL) %>% 
  ungroup() %>% 
  nest(-standtree)

BAI.fn <- function(x){gam(BAI ~ s(newring, k = 5), data = x)}
ringDens.fn <- function(x){gam(ringDens ~ s(newring, k = 5), data = x)}
LW.prop.fn <- function(x){gam(LW.prop ~ s(newring, k = 5), data = x)}
EW.dens.fn <- function(x){gam(EW.dens ~ s(newring, k = 5), data = x)}
LW.dens.fn <- function(x){gam(LW.dens ~ s(newring, k = 5), data = x)}
maxDens.fn <- function(x){gam(maxDens ~ s(newring, k = 5), data = x)}

df4 <- df3 %>% 
  mutate(BAImodel = map(data, BAI.fn),
         ringDensmodel = map(data, ringDens.fn),
         LW.propmodel = map(data, LW.prop.fn),
         EW.densmodel = map(data, EW.dens.fn),
         LW.densmodel = map(data, LW.dens.fn),
         maxDensmodel = map(data, maxDens.fn))

df5 <- df4 %>% 
  mutate(BAIPred = map2(BAImodel, data, predict),
         ringDensPred = map2(ringDensmodel, data, predict),
         LW.propPred = map2(LW.propmodel, data, predict),
         EW.densPred = map2(EW.densmodel, data, predict),
         LW.densPred = map2(LW.densmodel, data, predict),
         maxDensPred = map2(maxDensmodel, data, predict))

df6 <- df5 %>%
  unnest(data, BAIPred, ringDensPred, 
         LW.propPred, EW.densPred, LW.densPred, maxDensPred) %>% 
  filter(year >= 1962)

# This plot shows that geom_line now produces smoothed responses. This is a way to look at individual trees.
ggplot(filter(df6, PlotID == 7)) +
  geom_line(aes(x = year, y = ringDensPred, color = factor(Sample.HT))) +
  facet_wrap(~TreeID) +
  ylim(300, 860)

ggplot(df6) +
  stat_summary(aes(x = year, y = ringDens), fun.y = mean, geom = "line") +
  geom_smooth(aes(x = year, y = ringDens), method = "gam", formula = ~ s(k = 5)) +
  facet_wrap(~GSL) +
  ylim(300, 860)

# Graphs to show possible GSL effects on some responses at breast height
Fig5a <- ggplot(filter(df6, Sample.HT == 4.5)) +
  stat_summary(aes(x = year, y = BAIPred, color = factor(GSL)),
               fun.y = mean,
               geom = "line",
               size = 1) +
  ylab("BAI") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 16),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 24, margin = margin(r = 15))) +
  theme(legend.position = "none") 

Fig5b <- ggplot(filter(df6, Sample.HT == 4.5)) +
  stat_summary(aes(x = year, y = ringDensPred, color = factor(GSL)),
               fun.y = mean,
               geom = "line",
               size = 1) +
  ylab("RD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 16),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 24, margin = margin(r = 15))) +
  theme(legend.position = "none") 

Fig5c <- ggplot(filter(df6, Sample.HT == 4.5)) +
  stat_summary(aes(x = year, y = LW.propPred, color = factor(GSL)),
               fun.y = mean,
               geom = "line",
               size = 1) +
  ylab("LWP") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 24, margin = margin(r = 15))) +
  theme(legend.position = "none") 

Fig5aGrob <- ggplotGrob(Fig5a)
Fig5bGrob <- ggplotGrob(Fig5b)
Fig5cGrob <- ggplotGrob(Fig5c)
Fig5Grob <- rbind(Fig5aGrob, Fig5bGrob, Fig5cGrob, size = "first")
Fig5 <- grid.arrange(Fig5Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/RQ2.tiff", plot = Fig5,
       device = "tiff", 
       width = 16, height = 10, units = "cm")

Fig5Legend <- ggplot(df6) +
  geom_line(aes(x = year, y = ringDensPred, color = factor(GSL)), size = 1) +
  theme_bw() +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 24)) +
  labs(color = "GSL")

# Models using predicted values from GAMS would be very similar to a model just predicting the measured response, except climate and other high-freq variation could have an effect on the latter.

library(emmeans)

climWide <- read.csv("AnalysisReadyData/climWide_TWPaper.csv", 
                     header = TRUE) %>% 
  select(grow.yr, PDSIavg)
dfnew <- df2 %>% 
  left_join(climWide, by = c("year" = "grow.yr")) %>% 
  filter(year >= 1962)

# BAI
RQ2.1 <- lme(
  BAI ~ year + Sample.HT * factor(GSL),
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df2
)
furnival(RQ2.1)
Anova(RQ2.1, type = 3)
summary(RQ2.1)
df7$yhat <- predict(RQ2.1, level = 0)
bias(df7$yhat, df7$BAI)
plot(RQ2.1)

emmeans(RQ2.1, pairwise ~ GSL | Sample.HT)

# LWP
RQ2.2 <- lme(
  LW.prop ~ year + Sample.HT * factor(GSL),
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df2
)
furnival(RQ2.2)
Anova(RQ2.2, type = 3)
summary(RQ2.2)
df7$yhat <- predict(RQ2.2, level = 0)
bias(df7$yhat, df7$LW.propPred)
plot(RQ2.2)

emmeans(RQ2.2, pairwise ~ GSL | Sample.HT)

# RD
RQ2.3 <- lme(
  ringDens ~ year + Sample.HT * factor(GSL),
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df2
)
furnival(RQ2.3)
Anova(RQ2.3, type = 3)
summary(RQ2.3)
df7$yhat <- predict(RQ2.3, level = 0)
bias(df7$yhat, df7$ringDensPred)
plot(RQ2.3)

emmeans(RQ2.3, pairwise ~ GSL | Sample.HT)

# EWD
RQ2.4 <- lme(
  EW.densPred ~ year + Sample.HT * GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df6
)
furnival(RQ2.4)
Anova(RQ2.4, type = 3)
summary(RQ2.4)
df7$yhat <- predict(RQ2.4, level = 0)
bias(df7$yhat, df7$EW.densPred)
plot(RQ2.4)

# MXD
RQ2.5 <- lme(
  maxDensPred ~ year + Sample.HT * GSL,
  random = ~ 1 | TreeID,
  na.action = na.omit,
  data = df6
)
furnival(RQ2.5)
Anova(RQ2.5, type = 3)
summary(RQ2.5)
df7$yhat <- predict(RQ2.5, level = 0)
bias(df7$yhat, df7$maxDensPred)
