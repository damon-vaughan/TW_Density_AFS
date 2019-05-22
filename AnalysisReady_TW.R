# This script produces the files used in analysis in the TW project related to field measurements. Field taper measurements and acoustic velocity measurements

library(dplyr)
library(tidyr)

#############################
#Taylor Woods Taper
#############################

# Read in and create unique tree ID
twMeas <- read.csv2("IntermediateData/twMeas_2017.csv", header = TRUE,
                      colClasses = c("Date", rep("numeric", 23), "character")) %>%
  mutate(Tree.Key = paste(UnitID, TreeID, sep = "-"))
  
#######################################
# Taper Database
#######################################

twTaperData <- twMeas %>%
  select(-Date, -Notes, -CBH, -HM200_full, -HM200_16) 

# Gather wide to long, and fix the height column by removing X's
twTaperData <- twTaperData %>%
  gather(key = "Height", value = "Diameter", X0.5:X68) %>%
  mutate(Height = sub("X", "", Height))  %>%
  arrange(UnitID,TreeID)

# Change height to numeric- might want this as factor at some point
twTaperData$Height <- as.numeric(twTaperData$Height)

# Drop the NA's- might want these back at some point?
twTaperData <- twTaperData %>%
  na.omit()

# Write CSV
write.csv2(twTaperData,file = "AnalysisReadyData/twTaperData_2017.csv",
           row.names = FALSE)

##################################
# Acoustic Database
##################################

# remove columns that are not needed or provide tree-level info
twAcousticData <- twMeas %>%
  select(UnitID, TreeID, HM200_full, HM200_16, Tree.Key) 

# Drop the NA's- might want these back at some point?
twAcousticData <- twAcousticData %>%
  na.omit()

# Write CSV
write.csv2(twAcousticData,file = "AnalysisReadyData/twAcousticData_2017.csv",
           row.names = FALSE)