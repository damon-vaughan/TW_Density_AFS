library(readxl)
library(tidyverse)

# Individual tree data, from inventory list ----------------------------------------------------

TW <- read_excel("RawData/TW_2012_Marked_Data_newJune2016.xlsx",
                 sheet = "Marked trees",
                 col_names = c("PlotID", "TreeID", "DBH", "HT", "CBH", 
                               "HT.2016", "DBH.2016", "Velocity"),
                 skip = 2) %>%
  data.frame() %>% 
  mutate(TreeID = paste(PlotID, TreeID, sep = "-"))
str(TW)

# Drop non-sample trees. First create master tree list from Final WQ2017 Excel file
TreeList <- read_excel("RawData/TaylorWoods_Database_Summer17Final.xlsx") %>% 
  mutate(TreeID = paste(Unit, Tree, sep = "-")) %>% 
  select(TreeID)

TW.filtered <- left_join(TreeList,TW, by = "TreeID") %>%
  rename(Tree.Vel = Velocity)
# TW-13-176 did not get measured pre-felling. This is the only one.

# This plotsum block only used to help correct the missing tree; NOT used in actual plot summary data
plotsum <- TW.filtered %>% 
  group_by(PlotID) %>% 
  summarise(DBH = round(mean(DBH, na.rm = TRUE), 1), 
            HT = mean(HT, na.rm = TRUE), 
            CBH = round(mean(CBH, na.rm = TRUE), 1), 
            HT.2016 = round(mean(HT.2016, na.rm = TRUE), 1), 
            DBH.2016 = mean(DBH.2016, na.rm = TRUE))

TW.filtered[40,] <- c("13-176", 13, plotsum[13,2], plotsum[13,3], plotsum[13,4],
                   plotsum[13,5], plotsum[13,6], -999)

# Write to csv in Intermediate data folder
write.csv(TW.filtered, file = "IntermediateData/twTreeData.csv",
           row.names = FALSE)
tree <- read.csv("IntermediateData/twTreeData.csv",
                  header = TRUE)


# Plot-level data- uses other data sheet ----------------------------------
# As far as I can tell, not used in the paper

TW2 <- read.csv("RawData/TW_2015_NAU_data.csv", sep = ",",
                 col.names = c("SiteID", "Species", "PlotID", "Drop", 
                               "Plot.Num", "GSL", "TreeID", "DBH", "HT", 
                               "CBH", "Velocity", "Drop2", "Drop3", "Drop4",
                               "Drop5", "Drop6", "Drop7"),
                stringsAsFactors = FALSE) %>%
  filter(SiteID == "TW") %>%
  select(-SiteID, -Species, -Drop, -Drop2, -Drop3, -Drop4, -Drop5, 
         -Drop6, -Drop7) %>%
  data.frame() 

TW2$GSL <- as.numeric(TW2$GSL)

TW3 <- TW2 %>%
  arrange(PlotID, Plot.Num) %>%
  mutate(Tree.BA = DBH^2 * .005454,
         GSL = ifelse(GSL == 40, 30, GSL)) %>% 
  na.omit()

summary(TW3)

# Plot-level summaries
SubPlot.stats <- TW3 %>%
  group_by(PlotID, Plot.Num) %>%
  summarise(GSL = unique(GSL),
            BA = sum(Tree.BA) * 10,
            TPA = length(TreeID) * 10, 
            DBHmean = mean(DBH),
            DBHsd = sd(DBH, na.rm = T))
             
plotStats <- SubPlot.stats %>% 
  group_by(PlotID, GSL) %>% 
  summarise(BA = mean(BA),
            TPA = mean(TPA),
            DBH = mean(DBHmean),
            DBHsd = mean(DBHsd, na.rm = T)
            ) %>% 
  arrange(GSL)

plotStatsMetric <- plotStats %>% 
  mutate(GSL = round(GSL * 0.229568, 2),
         BA = round(BA * 0.229568, 2),
         TPH = round(TPA/0.404686, 1),
         DBH = round(DBH * 2.54, 2),
         DBHsd = round(DBHsd, 2)
         ) %>% 
  select(-TPA)

# write.csv(Plot.stats, "IntermediateData/TW_PlotStats.csv")
write.csv(plotStatsMetric, 
          "AnalysisReadyData/TWPaper_plotStats.csv", row.names = F)