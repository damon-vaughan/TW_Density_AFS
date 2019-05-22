library(readxl)
library(dplyr)

################################
#Taylor Woods Taper
###############################

#Read in and change column names. Generates warnings because of NA
twMeas <- read_excel("RawData/TaylorWoods_Database_Summer17Final.xlsx",
                     col_types = c("date", rep("numeric", 23), "text")) %>%
  data.frame()
colnames(twMeas) <- c("Date", "UnitID", "TreeID", "X0.5", "X4.5", "X8", "X12", 
                      "X16", "X20", "X24", "X28", "X32", "X36", "X40", "X44", 
                      "X48", "X52", "X56", "X60", "X64", "X68", "CBH", 
                      "HM200_full", "HM200_16", "Notes")

#Format columns
twMeas$Date <- as.Date(twMeas$Date)

#Output to CSV
write.csv2(twMeas,file = "IntermediateData/twMeas_2017.csv",
           row.names = FALSE)
