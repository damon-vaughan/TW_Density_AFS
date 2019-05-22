library(tidyverse)
library(data.table)

missing.rings <- fread("RawData/missing_rings.csv", sep = ",", header = TRUE)
maxyear <- read.csv("RawData/MaxYear.csv", header = TRUE, stringsAsFactors = FALSE)

filenames <- list.files("RawData/Scans_Spring2018/AnalyzedScans_DAT_SUM",
                        pattern = "dat$", full.names = TRUE)

df1 <- rbindlist(setNames(lapply(
  filenames, fread,  skip=40, fill = T, header = F, 
  col.names=c("junk","barktopith", "density","ring","EorL","volts")), 
  filenames),
  idcol = 'standtree') %>% 
  mutate(standtree = str_sub(standtree, 48, -5))

# Ring and latewood correction --------------------------------------------

df2 <- df1 %>% 
  left_join(maxyear, by = "standtree") %>% 
  nest(-standtree)

newBound2 <- function(x){
  df <- x %>%
    mutate(mm = max(barktopith) - barktopith + .01) %>%
    arrange(mm) %>% 
    mutate(ring = ring - (min(ring) - 1))
  
  df.1 <- df %>% 
    group_by(ring) %>% 
    summarise(start = min(mm), end = max(mm)) 
  
  df2 <- df %>% 
    left_join(df.1, by = c("ring")) %>% 
    group_by(ring) %>% 
    mutate(cutoff = min(density) + (max(density)-min(density)) * .8,
           latewood = ifelse(density > cutoff, 1, 0))
  
  df.2 <- df2 %>% 
    summarise(
      DistInRing1 = ifelse(max(latewood) == 1,
                           min(which(latewood == 1)) * .01 - .01,
                           0),
      DistInRing2 = ifelse(max(latewood) == 1,
                           max(which(latewood == 1)) * .01 - .01,
                           0),
      minmm = min(mm),
      boundary1 = DistInRing1 + minmm,
      boundary2 = DistInRing2 + minmm) 
  df.2[max(df.2$ring),"boundary2"] <- 1000
  
  df3 <- df2 %>% 
    ungroup() %>% 
    mutate(newring = as.numeric(cut(df2$mm, breaks = c(0,as.vector(
      df.2$boundary2))))) 
  
  df3 <- df3 %>% 
    left_join(df.2, by = c("newring" = "ring")) 
  
  df4 <- df3 %>% 
    ungroup() %>% 
    mutate(newEorL = ifelse(mm < boundary1, "EW", "LW"),
           year = unique(maxyear) - (max(newring) - newring)) %>% 
    select(mm, newring, density, cutoff, newEorL, year) %>% 
    as.data.frame()
  
  return(df4)
}

df3 <- df2 %>% 
  mutate(corrected = map(data, newBound2))

df4 <- df3 %>% 
  left_join(missing.rings, by = "standtree")

# Missing Rings -----------------------------------------------------------

missRing2002 <- function(x) {
  EW2002 <- x %>%
    group_by(year) %>%
    filter(year == 2003 & mm == min(mm)) %>%
    ungroup() %>%
    mutate(
      year = 2002,
      density = NA,
      cutoff = NA,
      newEorL = "EW"
    )
  LW2002 <- EW2002 %>%
    mutate(newEorL = "LW")
  x2 <- x %>%
    mutate(
      year = ifelse(year <= 2002, year - 1, year),
      newring = ifelse(year > 2002, newring + 1, newring)
    )
  x3 <- rbind(x2, EW2002, LW2002)
}

missRing1996 <- function(x) {
  EW1996 <- x %>%
    group_by(year) %>%
    filter(year == 1997 & mm == min(mm)) %>%
    ungroup() %>%
    mutate(
      year = 1996,
      density = NA,
      cutoff = NA,
      newEorL = "EW"
    )
  LW1996 <- EW1996 %>%
    mutate(newEorL = "LW")
  x2 <- x %>%
    mutate(
      year = ifelse(year <= 1996, year - 1, year),
      newring = ifelse(year > 1996, newring + 1, newring)
    )
  x3 <- rbind(x2, EW1996, LW1996)
}

missRing1956 <- function(x) {
  EW1956 <- x %>%
    group_by(year) %>%
    filter(year == 1957 & mm == min(mm)) %>%
    ungroup() %>%
    mutate(
      year = 1956,
      density = NA,
      cutoff = NA,
      newEorL = "EW"
    )
  LW1956 <- EW1956 %>%
    mutate(newEorL = "LW")
  x2 <- x %>%
    mutate(
      year = ifelse(year <= 1956, year - 1, year),
      newring = ifelse(year > 1956, newring + 1, newring)
    )
  x3 <- rbind(x2, EW1956, LW1956)
}

missRing1951 <- function(x) {
  EW1951 <- x %>%
    group_by(year) %>% 
    filter(year == 1952 & mm == min(mm)) %>% 
    ungroup() %>% 
    mutate(year = 1951, density = NA, cutoff = NA, newEorL = "EW")
  LW1951 <- EW1951 %>% 
    mutate(newEorL = "LW")
  x2 <- x %>% 
    mutate(year = ifelse(year <= 1951, year - 1, year),
           newring = ifelse(year > 1951, newring + 1, newring))
  x3 <- rbind(x2, EW1951, LW1951)
} 

df5 <- df4 %>% 
  mutate(fix2002 = ifelse(Y2002 == 1, map(corrected, missRing2002), 
                           corrected),
         fix2002 = map(fix2002, as.data.frame)) %>% 
  select(-data, -corrected) %>% 
  unnest(fix2002) %>% 
  select(-Y2002, -Y1996, -Y1956, -Y1951, -12) %>% 
  nest(-standtree) %>% 
  left_join(missing.rings, by = "standtree")

df6 <- df5 %>%
  mutate(fix1996 = ifelse(Y1996 == 1, map(data, missRing1996),
                          data),
         fix1996 = map(fix1996, as.data.frame)) %>%
  select(-data) %>%
  unnest(fix1996) %>%
  select(-Y2002,-Y1996,-Y1956,-Y1951,-12) %>%
  nest(-standtree) %>%
  left_join(missing.rings, by = "standtree")
  
df7 <- df6 %>%
  mutate(fix1956 = ifelse(Y1956 == 1, map(data, missRing1956),
                          data),
         fix1956 = map(fix1956, as.data.frame)) %>%
  select(-data) %>%
  unnest(fix1956) %>%
  select(-Y2002,-Y1996,-Y1956,-Y1951,-12) %>%
  nest(-standtree) %>%
  left_join(missing.rings, by = "standtree")

df8 <- df7 %>%
  mutate(fix1951 = ifelse(Y1951 == 1, map(data, missRing1951),
                          data),
         fix1951 = map(fix1951, as.data.frame)) %>%
  select(-data) %>%
  unnest(fix1951) %>%
  select(-Y2002,-Y1996,-Y1956,-Y1951,-12) 

# Summarise ---------------------------------------------------------------

FP <- df8 %>% 
  nest(-standtree)

sum_density <- function(x){
  # Summary data by tree and ring
  IR.dens <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(sg = mean(density, na.rm = T)) %>%  
    spread(newEorL, sg) %>% 
    rename(EWD = EW, LWD = LW)
  
  IR.width <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(width = ifelse(max(mm) != min(mm), 
                             max(mm) - min(mm) + .01, 0)) %>%  
    spread(newEorL, width) %>% 
    rename(EWW = EW, LWW = LW)
  
  distToMin <- x %>% 
    group_by(newring) %>% 
    na.omit %>% 
    summarise(distToMin = which.min(density) * .01 - .01)
  
  # mean density and ring width of each ring. Grouping sample, ring
  df1 <- x %>%
    group_by(newring) %>% 
    summarise(RW = ifelse(max(mm) != min(mm), max(mm) - min(mm) + .01, 0),
              RD = mean(density),
              minD = min(density),
              MXD = max(density),
              year = unique(year))
  
  # Joining sample, ring, and subring data. Separate ID into components
  df2 <- left_join(IR.dens, IR.width, by = "newring") %>% 
    left_join(distToMin, by = "newring") %>% 
    right_join(df1, by = "newring") 
  
  return(df2)
}

FP2 <- FP %>% 
  mutate(summarised = map(data, sum_density))

sumProfile <- FP2 %>% 
  unnest(summarised) %>% 
  separate(standtree, by = "-", 
           into = c("SiteID", "PlotID", "TreeID", "Position"), 
           remove = FALSE) %>% 
  unite("TreeID", PlotID, TreeID, sep = "-", remove = FALSE) %>% 
  select(standtree, SiteID, PlotID, TreeID, Position, newring, year,
         RW, EWW, LWW, distToMin, RD, minD, MXD, EWD, LWD)

write.csv(sumProfile, "IntermediateData/sumProfileInt.csv", 
          row.names = FALSE)

# Take summary profiles and produce dataframes for analysis ---------------

sumProfile <- read.csv("IntermediateData/sumProfileInt.csv", 
                       header = TRUE, stringsAsFactors = FALSE)

# Add mmRingStart
sumProfile2 <- sumProfile %>% 
  group_by(standtree) %>% 
  mutate(mmRingStart = cumsum(RW) - RW) %>% 
  ungroup()

# Position to height in meters
HT <- read.csv("IntermediateData/HT_to_join.csv")
sumProfile3 <- sumProfile2 %>% 
  left_join(HT, by = "Position") %>% 
  mutate(Sample.HT = Sample.HT * 0.3048)

# Add tree and plot level variables ---------------------------------------

# Tree level. For "13-176", no pre-felling data, so average others from the same plot
treedata <- read.csv("IntermediateData/twTreeData.csv", header = TRUE,
                     stringsAsFactors = FALSE) %>% 
  select(-Tree.Vel) 

sumProfile4 <- left_join(sumProfile3, treedata, by = "TreeID") %>% 
  select(-PlotID.y) %>% 
  rename(PlotID = PlotID.x)

# Plot level
plotdata <- read.csv("IntermediateData/TW_PlotStats.csv", header = TRUE) %>% 
  select(-X, -GSL)

sumProfile5 <- left_join(sumProfile4, plotdata, by = "PlotID")

# Treatment level
GSLdata <- read.csv("IntermediateData/GSL_to_join.csv")

sumProfile6 <- left_join(sumProfile5, GSLdata, by = "PlotID")

sumProfile7 <- sumProfile6 %>% 
  select(standtree, SiteID, PlotID, TreeID, Position, Sample.HT, 
         year, newring, RD, minD, MXD, EWD, LWD, RW, mmRingStart, 
         distToMin, EWW, LWW, DBH, HT, CBH, HT.2016, DBH.2016, BA, 
         QMD, AMD, TPA, GSL)

write.csv(sumProfile7, "IntermediateData/SumProfileFull.csv", 
          row.names = FALSE)

# Clean up the dataset ----------------------------------------------------

sumProfile <- read.csv("IntermediateData/SumProfileFull.csv", header = T, 
                       stringsAsFactors = F)

# Drop scans with poor correlations (From crossdating)
dropped <- read.csv("RawData/dropped_scans.csv", header = TRUE, 
                    stringsAsFactors = FALSE) 
dropped <- dropped$Scan

sumProfile2 <- sumProfile[!sumProfile$standtree %in% dropped,]

# The following sections dropped due to defects and curved rings
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

# Filter out first 10 rings and final ring, filter to heights <= 24 feet.

sumProfile8 <- sumProfile7 %>% 
  filter(Sample.HT <= 8) %>%
  group_by(standtree) %>% 
  filter(newring >= 10 & newring != max(newring)) %>% 
  ungroup() 

write.csv(sumProfile8, "AnalysisReadyData/sumProfile_AFS.csv", 
          row.names = FALSE)
