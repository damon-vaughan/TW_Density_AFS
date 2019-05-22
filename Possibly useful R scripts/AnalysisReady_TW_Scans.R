# Create analysis-ready data from hundreds of scans, that have already been fixed up by newBound in ScanCorrectionFull

library(tidyverse)
library(data.table)
source("missingRing.R")
source("sumDens.R")
source("EW_Correct_Fn.R")

# Create full dataframe from all scans... takes a while
filenames <- list.files("IntermediateData/CorrectedScans", full.names = TRUE)

d1 <- lapply(filenames, function(x){fread(x, header = TRUE)})

# Before summarizing, fix missing rings and max year
missing.rings <- read.csv("RawData/missing_rings.csv", header = TRUE, stringsAsFactors = FALSE)
maxyear <- read.csv("RawData/MaxYear.csv", header = TRUE, stringsAsFactors = FALSE)

# Apply missing ring corrections. missRingFull2 includes missing rings, giving them 0 widths and NA densities 
d2 <- lapply(d1, function(x){missRingFull2(x)})

fullProfile <- rbindlist(d2)
write.csv(fullProfile, "AnalysisReadyData/fullProfile.csv", row.names = FALSE)

# summarize each element of the list. sum_density4 works with missRingFull2 to give correct values to missing rings
d3 <- lapply(d2, function(x){sum_density4(x)})

sumProfile <- rbindlist(d3) %>% 
  separate(standtree, by = "-", 
           into = c("SiteID", "PlotID", "TreeID", "Position"), 
           remove = FALSE, convert = TRUE) %>% 
  unite("TreeID", PlotID, TreeID, sep = "-", remove = FALSE) %>% 
  select(standtree, SiteID, PlotID, TreeID, Position, newring, year,
         RW, EW.width, LW.width, distToMin, 
         ringDens, minDens, maxDens, EW.dens, LW.dens)

# Write to the Intermediate file- for Dendro analysis 
write.csv(sumProfile, file = "IntermediateData/sumProfileInt.csv",
          row.names = FALSE)

# Add variables and make corrections --------------------------------------
sumProfile <- read.csv("IntermediateData/sumProfileInt.csv", 
                       header = TRUE, stringsAsFactors = FALSE)
sumProfile <- read.csv("IntermediateData/sumProfileInt_purrr.csv", 
                             header = TRUE, stringsAsFactors = FALSE) 
sumProfile <- read.csv("IntermediateData/sumProfileInt_EFmethod.csv", 
                       header = TRUE, stringsAsFactors = FALSE) 

# Add mmRingStart
sumProfile2 <- sumProfile %>% 
  group_by(standtree) %>% 
  mutate(mmRingStart = cumsum(RW) - RW) %>% 
  ungroup()

# Position to height in feet
HT <- read.csv("IntermediateData/HT_to_join.csv")
sumProfile3 <- sumProfile2 %>% 
  left_join(HT, by = "Position") 

# If NA on EW width or LW width, both should be NA and intra-ring density should be NA. Likely results from fully decreasing density within the ring in a pith to bark direction. Could be compression wood or curved rings.
sumProfile4 <- sumProfile3 %>% 
  mutate(EW.width = ifelse(is.na(LW.width), NA, EW.width),
         LW.width = ifelse(is.na(EW.width), NA, LW.width),
         LW.dens = ifelse(is.na(EW.width), NA, LW.dens),
         EW.dens = ifelse(is.na(EW.width), NA, EW.dens))


# Add tree and plot level variables ---------------------------------------

# Tree level. For "13-176", no pre-felling data, so average others from the same plot
treedata <- read.csv("IntermediateData/twTreeData.csv", header = TRUE,
                     stringsAsFactors = FALSE) %>% 
  select(-Tree.Vel) 

plotsum <- treedata %>% 
  group_by(PlotID) %>% 
  summarise(DBH = round(mean(DBH, na.rm = TRUE), 1), 
            HT = mean(HT, na.rm = TRUE), 
            CBH = round(mean(CBH, na.rm = TRUE), 1), 
            HT.2016 = round(mean(HT.2016, na.rm = TRUE), 1), 
            DBH.2016 = mean(DBH.2016, na.rm = TRUE))
  
treedata[40,] <- c("13-176", 13, plotsum[13,2], plotsum[13,3], plotsum[13,4],
                   plotsum[13,5], plotsum[13,6])

write.csv(treedata, "AnalysisReadyData/treeData.csv", row.names = F)

sumProfile5 <- left_join(sumProfile4, treedata, by = "TreeID") %>% 
  select(-PlotID.y) %>% 
  rename(PlotID = PlotID.x)

# Plot level
plotdata <- read.csv("IntermediateData/TW_PlotStats.csv", header = TRUE) %>% 
  select(-X, -GSL)

sumProfile6 <- left_join(sumProfile5, plotdata, by = "PlotID")

# Treatment level
GSLdata <- read.csv("IntermediateData/GSL_to_join.csv")

sumProfile7 <- left_join(sumProfile6, GSLdata, by = "PlotID")

############################################
# Arrange and save
###########################################

sumProfile8 <- sumProfile7 %>% 
  select(standtree, SiteID, PlotID, TreeID, Position, Sample.HT, 
         year, newring, ringDens, minDens, maxDens,
         EW.dens, LW.dens, RW, mmRingStart, distToMin, EW.width, LW.width,
         DBH, HT, CBH, HT.2016, DBH.2016, BA, QMD, AMD, TPA, GSL)

write.csv(sumProfile8, "AnalysisReadyData/SumProfileFull_Jan.csv", 
          row.names = FALSE)
write.csv(sumProfile8, "AnalysisReadyData/SumProfileFull_purrr.csv", 
          row.names = FALSE)
write.csv(sumProfile8, "AnalysisReadyData/SumProfileFull_EFmethod.csv", 
          row.names = FALSE)
