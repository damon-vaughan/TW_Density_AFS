library(data.table)
library(tidyverse)
source("missingRing.R")

# Create full dataframe, apply missing ring correction, and save
filenames <- list.files("IntermediateData/CorrectedScans", full.names = TRUE)
d1 <- lapply(filenames, function(x){fread(x, header = TRUE)})

# Before summarizing, fix missing rings
missing.rings <- read.csv("RawData/missing_rings.csv", header = TRUE, stringsAsFactors = FALSE)
maxyear <- read.csv("RawData/MaxYear.csv", header = TRUE, stringsAsFactors = FALSE)
d2 <- lapply(d1, function(x){missRingFull(x)})

d3 <- d2 %>% rbindlist()

write.csv(d3, "IntermediateData/fullProfile.csv", row.names = FALSE)

fullProfile <- d3
fullProfile <- read.csv("IntermediateData/fullProfile.csv", sep = ",", 
                        header = TRUE)

test <- fullProfile %>% 
  filter(standtree == "TW-15-36-1") %>% 
  arrange(desc(mm)) %>% 
  mutate(year = 2015 - (max(newring) - newring))

test1 <- test[1:5000,]
test1a <- filter(test1, newEorL == "LW")
yearmm1 <- test1 %>% 
  group_by(newring) %>% 
  filter(mm == max(mm) & year %% 5 == 0)
# max mm is used because sometimes density is higher in the early part of ring
# label appears at the end of the year

test2 <- test[5000:10000,]
test2a <- filter(test2, newEorL == "LW")
yearmm2 <- test2 %>% 
  group_by(newring) %>% 
  filter(mm == max(mm) & year %% 5 == 0) 

test3 <- test[10000:15000,]
test3a <- filter(test3, newEorL == "LW")
yearmm3 <- test3 %>% 
  group_by(newring) %>% 
  filter(mm == max(mm) & year %% 5 == 0) 

test4 <- test[15000:20000,]
test4a <- filter(test4, newEorL == "LW")
yearmm4 <- test4 %>% 
  group_by(newring) %>% 
  filter(mm == max(mm) & year %% 5 == 0) 

test5 <- test[20000:25000,]
test5a <- filter(test5, newEorL == "LW")
yearmm5 <- test5 %>% 
  group_by(newring) %>% 
  filter(mm == max(mm) & year %% 5 == 0) 

ggplot(test1) +
  geom_line(aes(x = mm, y = density)) + 
  geom_point(aes(x = mm, y = density), data = test1a, color = "red") +
  scale_x_reverse() +
  geom_label(data = yearmm1, aes(x = mm, y = density, label = year),
             nudge_y = 100)

ggplot(test2) +
  geom_line(aes(x = mm, y = density)) +
  geom_point(aes(x = mm, y = density), data = test2a, color = "red") +
  scale_x_reverse() +
  geom_label(data = yearmm2, aes(x = mm, y = density, label = year),
             nudge_y = 100)

ggplot(test3) +
  geom_line(aes(x = mm, y = density)) +
  geom_point(aes(x = mm, y = density), data = test3a, color = "red") +
  scale_x_reverse() +
  geom_label(data = yearmm3, aes(x = mm, y = density, label = year),
             nudge_y = 50)

ggplot(test4) +
  geom_line(aes(x = mm, y = density)) +
  geom_point(aes(x = mm, y = density), data = test4a, color = "red") +
  scale_x_reverse() +
  geom_label(data = yearmm4, aes(x = mm, y = density, label = year),
            nudge_y = 50)

ggplot(test5) +
  geom_line(aes(x = mm, y = density)) +
  geom_point(aes(x = mm, y = density), data = test5a, color = "red") +
  scale_x_reverse() +
  geom_label(data = yearmm5, aes(x = mm, y = density, label = year),
             nudge_y = 50)
