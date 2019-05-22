# Add missing ring. Used in ScanCorrectionFull (summary modification part)

missRing <- function(x,y){
  df <- sumDens
  missing.year <- y
  newYear <- filter(df, year == missing.year) %>% 
    mutate(RW = 0, ringDens = 0, EW.dens = 0, LW.dens = 0, 
           EW.width = 0, LW.width = 0,
           newring = newring + 1)
  df2 <- df %>% 
    mutate(newring = ifelse(year > missing.year, newring+1, newring),
           year = ifelse(year <= missing.year, year-1, year)) %>% 
    rbind(newYear) %>% 
    arrange(newring)
}

# missRingFull is for full scans, not summaries. To be applied after newBound

missRingFull <- function(x){
  df <- x
  df2 <- left_join(df, maxyear, by = "standtree")
  df2 <- df2 %>% 
    mutate(year = unique(maxyear) - (max(newring) - newring))
  ringID <- semi_join(missing.rings, df2, by = "standtree")
  if(ringID$Y2002 == 1){
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 2002, year - 1, year),
             newring = ifelse(year > 2002, newring + 1, newring))} 
  if(ringID$Y1996 == 1){
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1996, year - 1, year),
             newring = ifelse(year > 1996, newring + 1, newring))} 
  if(ringID$Y1956 == 1){
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1956, year - 1, year),
             newring = ifelse(year > 1956, newring + 1, newring))} 
  if(ringID$Y1951 == 1){
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1951, year - 1, year),
             newring = ifelse(year > 1951, newring + 1, newring))} 
  if(ringID$Y1948 == 1){
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1948, year - 1, year),
             newring = ifelse(year > 1948, newring + 1, newring))} 
  return(df2)
  }

missRingFull2 <- function(x){
  df <- x
  df2 <- left_join(df, maxyear, by = "standtree")
  df2 <- df2 %>% 
    mutate(year = unique(maxyear) - (max(newring) - newring))
  ringID <- semi_join(missing.rings, df2, by = "standtree")
  if(ringID$Y2002 == 1){
    EW2002 <- df2 %>%
      group_by(year) %>% 
      filter(year == 2003 & mm == min(mm)) %>% 
      ungroup() %>% 
      mutate(year = 2002, density = NA, cutoff = NA, newEorL = "EW")
    LW2002 <- EW2002 %>% 
      mutate(newEorL = "LW")
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 2002, year - 1, year),
             newring = ifelse(year > 2002, newring + 1, newring))
    df2 <- rbind(df2, EW2002, LW2002)} 
  if(ringID$Y1996 == 1){
    EW1996 <- df2 %>%
      group_by(year) %>% 
      filter(year == 1997 & mm == min(mm)) %>% 
      ungroup() %>% 
      mutate(year = 1996, density = NA, cutoff = NA, newEorL = "EW")
    LW1996 <- EW1996 %>% 
      mutate(newEorL = "LW")
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1996, year - 1, year),
             newring = ifelse(year > 1996, newring + 1, newring))
    df2 <- rbind(df2, EW1996, LW1996)} 
  if(ringID$Y1956 == 1){
    EW1956 <- df2 %>%
      group_by(year) %>% 
      filter(year == 1957 & mm == min(mm)) %>% 
      ungroup() %>% 
      mutate(year = 1956, density = NA, cutoff = NA, newEorL = "EW")
    LW1956 <- EW1956 %>% 
      mutate(newEorL = "LW")
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1956, year - 1, year),
             newring = ifelse(year > 1956, newring + 1, newring))
    df2 <- rbind(df2, EW1956, LW1956)} 
  if(ringID$Y1951 == 1){
    EW1951 <- df2 %>%
      group_by(year) %>% 
      filter(year == 1952 & mm == min(mm)) %>% 
      ungroup() %>% 
      mutate(year = 1951, density = NA, cutoff = NA, newEorL = "EW")
    LW1951 <- EW1951 %>% 
      mutate(newEorL = "LW")
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1951, year - 1, year),
             newring = ifelse(year > 1951, newring + 1, newring))
    df2 <- rbind(df2, EW1951, LW1951)} 
  if(ringID$Y1948 == 1){
    EW1948 <- df2 %>%
      group_by(year) %>% 
      filter(year == 1949 & mm == min(mm)) %>% 
      ungroup() %>% 
      mutate(year = 1948, density = NA, cutoff = NA, newEorL = "EW")
    LW1948 <- EW1948 %>% 
      mutate(newEorL = "LW")
    df2 <- df2 %>% 
      mutate(year = ifelse(year <= 1948, year - 1, year),
             newring = ifelse(year > 1948, newring + 1, newring))
    df2 <- rbind(df2, EW1948, LW1948)} 
  return(df2)
}