# Takes a dataframe of a single density profile and creates a dataframe summarized by year

###############################
# sum_density is for individual scans, as in ScanCorrectionFull.R
###############################

sum_density <- function(x){
# Summary data by tree and ring
df1a <- corScan %>% 
  group_by(newring, newEorL) %>% 
  summarise(sg = mean(density, na.rm = T)) %>% 
  spread(newEorL, sg)

# mean density and ring width of each ring. Grouping sample, ring
df1b <- corScan %>%
  group_by(standtree, newring) %>% 
  summarise(RW = max(mm) - min(mm),
            ringDens = mean(density),
            minDens = min(density),
            maxDens = max(density)) 

# create year column. Grouping sample
df1c <- df1b %>% 
  group_by(standtree) %>% 
  mutate(year = 2017 - (max(newring) - newring))

# Joining sample, ring, and subring data. Separate ID into components
df2 <- left_join(df1a, df1c, by = "newring") %>%
  separate(standtree, c("SiteID", "PlotID", "TreeID", "Position"), 
           remove = FALSE)

# Finalize dataframe. Drop first ring, last ring, and ungroup
sumDens <- df2 %>% 
  group_by(standtree) %>% 
  # filter(newring != "1" & newring != max(newring)) %>% 
  na.omit() %>% 
  ungroup(sumDens)

sumDens <<- sumDens
}

############################
# what is sum_density2 for?
###########################

sum_density2 <- function(x){
  # Summary data by tree and ring
  df1a <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(sg = mean(density, na.rm = T)) %>% 
    spread(newEorL, sg)
  
  # mean density and ring width of each ring. Grouping sample, ring
  df1b <- x %>%
    group_by(standtree, newring) %>% 
    summarise(RW = max(mm) - min(mm),
              ringDens = mean(density),
              minDens = min(density),
              maxDens = max(density)) 
  
  # create year column. Grouping sample
  df1c <- df1b %>% 
    group_by(standtree) %>% 
    mutate(year = 2017 - (max(newring) - newring))
  
  # Joining sample, ring, and subring data. Separate ID into components
  df2 <- left_join(df1a, df1c) %>%
    separate(standtree, c("SiteID", "PlotID", "TreeID", "Position"), 
             remove = FALSE)
  
  # Finalize dataframe. Drop first ring, last ring, and ungroup
  sumDens <- df2 %>% 
    group_by(standtree) %>% 
    filter(newring != "1" & newring != max(newring)) %>% 
    na.omit() %>% 
    ungroup(sumDens)
  
  sumDens <<- sumDens
}

#########################################
# sum_density3 is for use on full dataframe. Also works on individual scans
#######################################

sum_density3 <- function(x){
  # Summary data by tree and ring
  IR.dens <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(sg = mean(density, na.rm = T)) %>%  
    spread(newEorL, sg) %>% 
    rename(EW.dens = EW, LW.dens = LW)
  
  IR.width <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(width = max(mm) - min(mm) + .01) %>%  
    spread(newEorL, width) %>% 
    rename(EW.width = EW, LW.width = LW)
  
  # mean density and ring width of each ring. Grouping sample, ring
  df1 <- x %>%
    group_by(standtree, newring) %>% 
    summarise(RW = max(mm) - min(mm) + .01,
              ringDens = mean(density),
              minDens = min(density),
              maxDens = max(density),
              year = unique(year))

  # Joining sample, ring, and subring data. Separate ID into components
  df2 <- left_join(IR.dens, IR.width, by = "newring") %>% 
    left_join(df1, by = "newring") %>% 
    separate(standtree, c("SiteID", "PlotID", "TreeID", "Position"), 
             remove = FALSE)
  
  # Finalize dataframe
  sumDens <- df2 %>% 
    group_by(standtree) %>% 
    ungroup(sumDens)
  
  sumDens <<- sumDens
}

sum_density4 <- function(x){
  # Summary data by tree and ring
  IR.dens <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(sg = mean(density, na.rm = T)) %>%  
    spread(newEorL, sg) %>% 
    rename(EW.dens = EW, LW.dens = LW)
  
  IR.width <- x %>% 
    group_by(newring, newEorL) %>% 
    summarise(width = ifelse(max(mm) != min(mm), max(mm) - min(mm) + .01, 0)) %>%  
    spread(newEorL, width) %>% 
    rename(EW.width = EW, LW.width = LW)
  
  distToMin <- x %>% 
    group_by(newring) %>% 
    na.omit %>% 
    summarise(distToMin = which.min(density) * .01 - .01)

  # mean density and ring width of each ring. Grouping sample, ring
  df1 <- x %>%
    group_by(standtree, newring) %>% 
    summarise(RW = ifelse(max(mm) != min(mm), max(mm) - min(mm) + .01, 0),
              ringDens = mean(density),
              minDens = min(density),
              maxDens = max(density),
              year = unique(year))
  
  # Joining sample, ring, and subring data. Separate ID into components
  df2 <- left_join(IR.dens, IR.width, by = "newring") %>% 
    left_join(distToMin, by = "newring") %>% 
    right_join(df1, by = "newring") %>% 
    separate(standtree, c("SiteID", "PlotID", "TreeID", "Position"), 
             remove = FALSE)
  
  # Finalize dataframe
  sumDens <- df2 %>% 
    group_by(standtree) %>% 
    ungroup(sumDens)
  
  sumDens <<- sumDens
}
