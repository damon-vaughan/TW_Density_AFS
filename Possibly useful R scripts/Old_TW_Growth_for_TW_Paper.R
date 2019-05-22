# Use full summary dataframe including first and last ring, but drop the scans that don't correlate.

sumProfile <- read.csv("AnalysisReadyData/sumProfileFull_Jan.csv", 
                       header = TRUE,
                       stringsAsFactors = FALSE)

dropped <- read.csv("RawData/dropped_scans.csv", header = TRUE, stringsAsFactors = FALSE) 
dropped <- dropped$Scan

sumProfile2 <- sumProfile[!sumProfile$standtree %in% dropped,]

# climWide <- read.csv("AnalysisReadyData/climWide.csv", header = TRUE)
# df2 <- left_join(sumProfile2, climWide, by = c("year" = "grow.yr"))

# This df has only breast height, and 2017 is excluded

df2 <- sumProfile2 %>% 
  mutate(LW.prop = LW.width/RW,
         BAI = pi * ((mmRingStart + RW)^2 - mmRingStart^2)/100,
         EWAI = pi * ((mmRingStart + EW.width)^2 - mmRingStart^2)/100,
         LWAI = pi * ((mmRingStart + RW)^2 - (mmRingStart + RW - LW.width)^2)/100,
         LWAI.prop = LWAI/BAI) %>% 
  filter(Sample.HT == 4.5 & year < 2017) %>% 
  na.omit

# Standard errors for geom_ribbon

df2.sum <- df2 %>% 
  group_by(GSL, year) %>% 
  summarise(mean = mean(BAI, na.rm = TRUE), sd = sd(BAI, na.rm = TRUE)) %>% 
  na.omit

ggplot(df2.sum) +
  geom_line(aes(x = year, y = mean, color = factor(GSL))) +
  geom_ribbon(aes(x = year, ymin = mean - sd, ymax = mean + sd, fill = factor(GSL)), alpha = .4)

# cumulative area graph

df2 <- df2 %>% 
  group_by(standtree) %>% 
  mutate(BAI = pi * ((mmRingStart + RW)^2 - mmRingStart^2)/100,
         CumulativeBA = cumsum(BAI))

ggplot(df2) +
  geom_smooth(aes(x = year, y = CumulativeBA, color = factor(GSL))) +
  facet_wrap(~Sample.HT) +
  geom_vline(xintercept = 1962)

ggplot(filter(df2, TreeID == "16-44", Sample.HT == 4.5)) +
  geom_point(aes(x = year, y = CumulativeBA)) +
  # facet_wrap(~TreeID) +
  geom_vline(xintercept = 1962)

ggplot(filter(df2, Sample.HT == 4.5 & TreeID == "8-657")) +
  geom_point(aes(x = year, y = CumulativeBA)) +
  # facet_wrap(~TreeID) +
  geom_vline(xintercept = 1962)

ggplot(filter(df2, Sample.HT == 4.5)) +
  geom_smooth(aes(x = year, y = CumulativeBA, color = factor(GSL))) +
  # facet_wrap(~GSL) +
  geom_vline(xintercept = 1962)

# Check to see growth differences prior to 1962

pre62 <- df2 %>% 
  filter(year <= 1962) %>% 
  filter(Sample.HT == 4.5)

pre62mod <- lm(BAI ~ PDSIavg, data = pre62)
summary(pre62mod)
anova(mod1962)

# modeling growth response to thinning- cumulative BA

thin.yr <- 1962

pre.post.thin <- df2 %>% 
  filter(year >= thin.yr - 5 & year <= thin.yr + 15 & Sample.HT == 4.5) %>% 
  mutate(pre.post = ifelse(year <= thin.yr, "pre", "post")) %>%
  ungroup() %>% 
  select(TreeID, year, GSL, PlotID, BAI, LW.prop)

mod62 <- lm(BAI ~ pre.post, data = pre.post.thin)
summary(mod62)

ggplot(pre.post.thin) +
  stat_summary(aes(x = year, y = LW.prop, color = factor(GSL)), fun.y = mean, geom = "line")

# piecewise linear
pre.post.thin <- pre.post.thin %>% 
  mutate(knot.ind = ifelse(pre.post.thin$year > thin.yr, (pre.post.thin$year - thin.yr), 0))

mod <- lmList(CumulativeBA ~ year + knot.ind | TreeID, data = pre.post.thin)

TreeID <- unique(pre.post.thin$TreeID)

coefs <- coef(mod)[,2:3] %>% 
  data.frame() %>% 
  rename(slope.pre = year, slope.change = knot.ind)
pvals <- summary(mod)$coefficients[,4,2:3] %>% 
  data.frame() %>% 
  rename(slope.pre.p = year, slope.change.p = knot.ind)

results <- cbind(TreeID, coefs, pvals) %>% 
  mutate(slope.pre = ifelse(slope.pre.p < 0.1, slope.pre, 0),
         slope.change = ifelse(slope.change.p < 0.1, slope.change, 0))

TreeGSL <- pre.post.thin %>% 
  select(TreeID, PlotID, GSL) %>% 
  unique()

results <- results %>% 
  left_join(TreeGSL, by = "TreeID")

# Change in the "post" period. Positive indicates a release, negative means crowding
ggplot(results) +
  geom_bar(aes(x = TreeID, y = slope.change, fill = factor(GSL)), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = .25))

resultsum <- results %>% 
  filter(slope.change != 0) %>% 
  group_by(PlotID) %>% 
  summarise(slope.change = mean(slope.change), GSL = unique(GSL))

ggplot(resultsum) +
  geom_bar(aes(x = PlotID, y = slope.change, fill = factor(GSL)), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = .25))  

resultsum2 <- results %>% 
  filter(slope.change != 0) %>% 
  group_by(GSL) %>% 
  summarise(slope.change = mean(slope.change))

ggplot(resultsum2) +
  geom_bar(aes(x = GSL, y = slope.change, fill = factor(GSL)), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = .25))  

# Values from pre-post. This is in absolute terms, so it cannot be negative
results2 <- results %>% 
  gather(key = "pre.post", value = "value", c("year", "slope.post"))

ggplot(results2) +
  geom_bar(aes(x = TreeID, y = value, color = pre.post), stat = "identity", position = "dodge")

# Modeling growth response to climate - Annual BAI 

thin.yr <- 2002

pre.post.thin <- df2 %>% 
  filter(year >= thin.yr - 5 & year <= thin.yr + 7 & Sample.HT == 4.5) %>% 
  mutate(pre.post = ifelse(year <= thin.yr, "pre", "thinned"),
         years.after = ifelse(year > thin.yr, year - thin.yr, 0)) %>%
  ungroup() %>% 
  select(BAI, TreeID, year, PDSIavg, prcptot, RW, pre.post, years.after, GSL)

# mod <- lmList(BAI ~ PDSIavg + years.after | factor(TreeID), data = pre.post.thin)
# PDSIeffect <- summary(mod)$coefficients[,c(1,4),2]

mod <- lm(BAI ~ PDSIavg + factor(pre.post) + factor(GSL), data = pre.post.thin)
anova(mod)
summary(mod)

mm1 <- lmer(BAI ~ PDSIavg + factor(pre.post) + (1|GSL), data = pre.post.thin)
summary(mm1)
anova(mm1)

ggplot(filter(pre.post.thin, TreeID == "9-317")) +
  geom_point(aes(x = year, y = BAI))

ggplot(filter(pre.post.thin, GSL == 150)) +
  geom_point(aes(x = year, y = BAI)) +
  geom_vline(xintercept = thin.yr) +
  facet_wrap(~TreeID)

ggplot(df2) +
  geom_point(aes(x = year, y = BAI)) 

# moving to full dataframe

fullat4.5 <- df2 %>% 
  filter(Sample.HT == 4.5) %>% 
  ungroup() 

ggplot(filter(fullat4.5, GSL == 150)) +
  geom_point(aes(x = year, y = BAI)) +
  facet_wrap(~TreeID) +
  geom_vline(xintercept = 1962)

# post-1962
post62 <- df2 %>% 
  filter(year >= 1962) %>% 
  filter(Sample.HT == 4.5)

mod <- lm(RW ~ PDSIavg * factor(GSL) + mmRingStart, data = post62)
anova(mod)
summary(mod)

# mod <- lm(LW.width/RW ~ PDSIavg + factor(GSL) + mmRingStart, data = post62)
# summary(mod)

ggplot(filter(post62, GSL == 60)) +
  geom_smooth(aes(x = year, y = BAI)) +
  facet_wrap(~TreeID) 

ggplot(filter(post62)) +
  geom_smooth(aes(x = PDSIavg, y = RW, color = factor(GSL))) 

ggplot(filter(post62)) +
  geom_smooth(aes(x = year, y = RW, color = factor(GSL))) 

post62 <- df2 %>% 
  filter(year >= 1962) 

mod <- lm(BAI ~ PDSIavg + factor(GSL) + Sample.HT, data = post62)
summary(mod)


# TW RQ #2 (or #3?)- Do ring width increases also affect latewood proportion? --------------------------

df3 <- df2 %>% 
  filter(GSL == 30 | GSL == 60)

ggplot(df3) +
  geom_smooth(aes(x = year, y = BAI, color = factor(GSL)))

ggplot(df2) +
  stat_summary(aes(x = year, y = RW, color = factor(GSL)), fun.y = mean, geom = "line")

df4 <- df2 %>% 
  filter(newring >= 20) 

ggplot(df4) +
  stat_summary(aes(x = year, y = LW.dens, color = factor(GSL)), fun.y = mean, geom = "line")
# Latewood density goes up after thinning, but doesn't separate by treatment. 
ggplot(df4) +
  stat_summary(aes(x = year, y = EW.dens, color = factor(GSL)), fun.y = mean, geom = "line")
# Earlywood density spikes in drought years. Supports the slide I showed in New Frontiers. Result of less cell expansion, possibly false rings but with my method false ring would trigger latewood classification?
# No effect of 1962 thinning on earlywood density
