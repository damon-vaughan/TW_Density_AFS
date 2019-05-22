library(tidyverse)
library(dplR)

# Try funciton "column_to_rownames" in dplyr

sumProfile <- read.csv("IntermediateData/sumProfileInt.csv", 
                       sep = ",", stringsAsFactors = FALSE,
                       header = TRUE)

dropped <- read.csv("RawData/dropped_scans.csv", header = TRUE, stringsAsFactors = FALSE) 
dropped <- dropped$Scan

sumProfile2 <- sumProfile[!sumProfile$standtree %in% dropped,]

d <- sumProfile2 %>%
  filter(sumProfile3$PlotID == 3) %>%
  select(TreeID,year,Position,RW) %>%
  unite(NewID, TreeID, sep = "-", Position, remove = TRUE) %>%
  spread(key = NewID, value = RW) 

rownames(d) <- d$year

d2 <- d %>% 
  select(-year)

# Raw data as rwl file is available from https://www.ncdc.noaa.gov/paleo-search/?dataTypeId=18
corr.20 <- corr.rwl.seg(d2, seg.length = 20, bin.floor = 0, pcrit = .1)

which(corr.20$overall[,1] <= 0.35)
mean(corr.20$overall[,1])
# Mean inter-series correlation is .6597


# Unchanged below ---------------------------------------------------------

df <- d2

flagged <- df$`1_188_4`
names(flagged) <- rownames(df)
df$`1_188_4` <- NULL

corr.series.seg(d2, series = flagged, seg.length = 10, 
                bin.floor = 0, pcrit = .1)

ccf.series.rwl(rwl = df, series = flagged, seg.length = 10, bin.floor = 0)

ggplot(sumProfile[sumProfile$standtree == "TW-1-146-45",]) +
  geom_bar(aes(x = year, y = RW), stat = "identity") +
  geom_vline(xintercept = 2002) +
  geom_vline(xintercept = 1996) +
  geom_vline(xintercept = 1989) +
  geom_vline(xintercept = 1977) 
+
  geom_vline(xintercept = 2012)

