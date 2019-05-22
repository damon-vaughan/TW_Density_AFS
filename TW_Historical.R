# Data from Kristen- this includes pre-thinning data measured on all trees, not sub-plots, up until 2003. Missing post-thinning data, pre-thin 1962, and anything after 2003. I've checked this against Ronco reported basal areas and its almost exactly the same.

# Also available here are individual tree DBH (to make average), height, and crown length. Maybe could use this in the future but not included right now.

# an "a" after year indicates pre-thin, "b" indicates post-thin

library(tidyverse)
library(readxl)

TW_1962 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1962", col_names = TRUE) %>% 
  mutate(year = "Y1962b")
TW_1967 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1967", col_names = TRUE) %>% 
  mutate(year = "Y1967")
TW_1972 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1972", col_names = TRUE) %>% 
  mutate(year = "Y1972a")
TW_1977 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1977", col_names = TRUE) %>% 
  mutate(year = "Y1977") 
TW_1982 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1982", col_names = TRUE) %>% 
  mutate(year = "Y1982a")
TW_1987 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1987", col_names = TRUE) %>% 
  mutate(year = "Y1987")
TW_1992 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1992", col_names = TRUE) %>% 
  mutate(year = "Y1992a")
TW_1998 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "1998", col_names = TRUE) %>% 
  select(-Burn) %>% 
  mutate(year = "Y1998")
TW_2003 <- read_excel("RawData/Taylor Woods Master.xls", sheet = "2003", col_names = TRUE) %>% 
  mutate(year = "Y2003a") %>% 
  mutate(Diam. = ifelse(Plot == 8 & Tree == 543, Diam./10, Diam.))

TW_History_Long <- rbind(TW_1962, TW_1967, TW_1972, TW_1977, TW_1982, TW_1987, TW_1992, TW_1998, TW_2003) %>% 
  select(year, Plot, GSL, Tree, Diam.) %>% 
  mutate(Diam. =  Diam./10,
         BA = Diam.^2 * .005454) %>% 
  filter(Tree != 543 | Plot != 8 | year != "Y2003a")

TW_Unit_Acreage <- data.frame(Plot = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), 
                              Acres = c(.8,.8,.8,1,.8,.82,1.24,1,.75,.75,.75,.75,.75,.75,.75,.75,.75,1))

TW_Sum_Long <- TW_History_Long %>% 
  group_by(year, Plot) %>% 
  summarise(TotBA = sum(BA, na.rm = TRUE)) %>% 
  na.omit() %>% 
  left_join(TW_Unit_Acreage, by = "Plot") %>% 
  mutate(BA.acre = TotBA/Acres) %>% 
  select(-TotBA, -Acres)

TW_History <- TW_Sum_Long %>% 
  spread(key = year, value = BA.acre) 

check <- TW_History_Long %>% 
  filter(year == "Y2003a") %>%
  na.omit() %>% 
  group_by(GSL, Plot) %>% 
  summarise(meanDia = mean(Diam., na.rm = T),
            count = length(Diam.),
            maxD = max(Diam.))
  
# Supplement above df with data from Ronco. Addings 1962a, 1972b, 1982b

RoncoData <- read.csv("RawData/RoncoData.csv", header = TRUE, stringsAsFactors = FALSE)

TW_History2 <- TW_History %>% 
  left_join(RoncoData, by = c("Plot" = "PlotID")) %>% 
  select(Plot, GSL, Y1962a, Y1962b, Y1967, Y1972a, Y1972b, Y1977, Y1982a, Y1982b, Y1987, Y1992a, Y1998, Y2003a)

# Add Dave's 2015 measurements

TW_2015 <- read.csv("RawData/TW_2015_NAU_data.csv", header = TRUE, stringsAsFactors = FALSE) %>% 
  filter(site == "TW") %>%
  select(Plot = stand, subplot = plot, treatment, tree, dbh_in) %>% 
  mutate(BA = dbh_in^2 * .005454)

SubPlot.sum <- TW_2015 %>%
  group_by(Plot, subplot) %>%
  summarise(BA = sum(BA) * 10) %>% 
  filter(Plot != 19 & Plot != 22)  

Plot.sum <- SubPlot.sum %>% 
  group_by(Plot) %>% 
  summarise(Y2015a = mean(BA))

TW_History3 <- TW_History2 %>% 
  left_join(Plot.sum, by = "Plot")

# Add post-thin basal areas where they are not measured (1992, 2003, 2017). These are assumed based on average pre-thinning tree diameter. Should be based on avg post thin-diameter, but this data not available.

TW_AvgDia <- rbind(TW_1962, TW_1967, TW_1972, TW_1977, TW_1982, TW_1987, TW_1992, TW_1998, TW_2003) %>% 
  filter((year == "Y1992a" | year == "Y2003a") & Diam. > 0) %>% 
  select(year, Plot, Diam.) %>% 
  mutate(Diam. = Diam./10) %>% 
  group_by(year, Plot) %>% 
  summarise(Diam. = mean(Diam., na.rm = TRUE)) %>% 
  na.omit() %>% 
  ungroup()

TW_2015_Dia <- TW_2015 %>% 
  group_by(Plot) %>% 
  summarise(Diam. = mean(dbh_in)) %>% 
  filter(Plot != 19 & Plot != 22) %>% 
  mutate(year = "Y2015a")

GSLcol <- RoncoData %>% 
  select(GSL, PlotID)

# Use these average diameters, along with table from Myers 1967, to determine post-thinning density. Almost all plots average over 10" so it is rarely different from the GSL
TW_AvgDia2 <- rbind(TW_AvgDia, TW_2015_Dia) %>% 
  left_join(GSLcol, by = c("Plot" = "PlotID"))

TW_PostThins <- RoncoData %>% 
  select(GSL, PlotID) %>% 
  mutate(Y1992b = GSL, Y2003b = GSL, Y2017b = GSL) %>% 
  mutate(Y1992b = ifelse(PlotID == 4, 141.4, 
                         ifelse(PlotID == 5, 117,
                                ifelse(PlotID == 8, 141.4,
                                       ifelse(PlotID == 10, 141.4,
                                       Y1992b))))) %>% 
  select(-GSL)

TW_History4 <- left_join(TW_History3, TW_PostThins, by = c("Plot" = "PlotID")) %>% 
  select(Plot:Y1992a, Y1992b, Y1998, Y2003a, Y2003b, Y2015a, Y2017b) 

TW_History_Metric <- round(TW_History4 * 0.229568, 1)

write.csv(TW_History_Metric, "IntermediateData/TW_History_BA.csv", row.names = FALSE)
TW_History <- read.csv("IntermediateData/TW_History_BA.csv", header = TRUE)

# Make long
TW_History <- TW_History_Metric

TW_Long <- TW_History %>% 
  gather("YearID", "BA", Y1962a:Y2017b) 

TW_Long2 <- TW_Long %>% 
  separate(YearID, c(1,5), into = c("Y", "year", "pre.post"), convert = TRUE) %>% 
  select(-Y) %>% 
  mutate(year = ifelse(pre.post == "b", year + .5, year))

str(TW_Long2)

write.csv(TW_Long2, "AnalysisReadyData/TW_Growth_ForPaper.csv", row.names = FALSE)

# Plot of BA over time for the 6 treatments
ggplot(TW_Long2) +
  # geom_line(aes(x = year, y = BA, color = factor(Plot))) +
  stat_summary(aes(x = year, y = BA, color = factor(GSL)), fun.y = mean, geom = "line", size = 1) +
  geom_hline(yintercept = c(6.9, 13.8, 18.4, 23, 27.6, 34.4), color = c("orangered2", "gold3", "springgreen3", 
                                                                 "turquoise3", "dodgerblue2", "magenta")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab(expression(paste("Basal Area ", "(ft" ^ "2 ","ac" ^ "-1", ")")))

startval <- TW_History %>%
  group_by(GSL) %>% 
  summarise(mean = mean(Y1962a))
mean(startval$mean) * 0.229568

# Percent cut
TW_Cut <- TW_History %>% 
  mutate(cut1962 = 1 - Y1962b/Y1962a,
         cut1972 = 1 - Y1972b/Y1972a,
         cut1982 = 1 - Y1982b/Y1982a,
         cut1992 = 1 - Y1992b/Y1992a,
         cut2003 = 1 - Y2003b/Y2003a,
         cut2017 = 1 - Y2017b/Y2015a) %>% 
  select(Plot, GSL, cut1962, cut1972, cut1982, cut1992, cut2003, cut2017)

TW_Cut_Long <- TW_Cut %>% 
  gather("cutID", "percent.cut", cut1962:cut2017)
