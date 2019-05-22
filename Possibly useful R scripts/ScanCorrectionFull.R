# Adjusts QTR ring boundaries with "newBound", summarizes by year with "sum_density", adds missing rings with "missRing", and adds or replaces values in full summary dataframe.
# Based on "Cleanup_TW_Scans" which is now in the archive

library(tidyverse)
library(data.table)
source("NewRing_NewBound.R")
source("sumDens.R")
source("missingRing.R")

rawScan <- fread("RawData/Scans_Spring2018/AnalyzedScans_DAT_SUM/TW-17-19-6.dat",
            skip=40, 
            fill = T, 
            header = F, 
            col.names=c("junk","barktopith", "density",
                        "ring","EorL","volts")) %>% 
  na.omit() %>% 
  mutate(standtree = "TW-17-19-6") %>% 
  select(standtree, c(1:6))

missing.rings <- fread("RawData/missing_rings.csv", sep = ",", header = TRUE)
maxyear <- read.csv("RawData/MaxYear.csv", header = TRUE, stringsAsFactors = FALSE)

# Fix LW/ring boundaries- sourced from "NewRing_NewBound"
newBound(rawScan)

write.table(corScan, file = 'IntermediateData/CorrectedScans/TW-17-19-6', 
            row.names = F, col.names = T, sep = "\t")

###################
# Apply fixes to the already aggregated summary table
###################

corScan <- missRingFull2(corScan)
sum_density4(corScan)

# Use below line if the series does not start at 2017. This corrects AFTER the sum_density, which already dropped 2017. So the subtraction part might not make sense at first

# sumDens <- sumDens %>%
  # mutate(year = year - 3)

# Adds missing rings. Calls function and relies on manual entry into textfile. This was so hard to do!!
# ringID <- semi_join(missing.rings,sumDens, by = "standtree")

# if(ringID$Y2002 == 1){
#   sumDens <- missRing(sumDens, 2002)} 
# 
# if(ringID$Y1996 == 1){
#   sumDens <- missRing(sumDens, 1996)} 
# 
# sumDens2 <- sumDens %>%
#   left_join(missing.rings, by = "standtree")
  
sumDens2 <- sumDens %>% 
  select(newring,EW.dens,LW.dens,EW.width,LW.width,standtree,SiteID,PlotID,
         TreeID,Position,RW,ringDens,minDens,maxDens,year)

# Add new scan to dataframe and replace old
sumProfile <- read.csv("IntermediateData/sumProfileInt.csv", sep = ",",
                       header = TRUE,
                       stringsAsFactors = FALSE)

sumProfile2 <- sumProfile %>%
  filter(standtree != unique(corScan$standtree)) %>% 
  rbind(sumDens2) %>% 
  arrange(standtree)

write.csv(sumProfile2,
          file = "IntermediateData/sumProfileInt.csv", 
          row.names = FALSE)
