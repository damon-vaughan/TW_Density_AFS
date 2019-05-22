# Takes roughly delineated scans from QMS and:
# 1. Arranges by mm, bark to pith
# 2. Identify EW/LW boundary as 80% of the max-min
# 3. Identify new ring boundary as the end of the previous ring's LW
# 4. Join to original dataframe

newBound <- function(x){
  # Step 1: Arrange by mm, delete non-wood measurements, flips readings so they    are pith to bark. 
  
  df <- rawScan %>%
    mutate(mm = max(barktopith) - barktopith + .01) %>%
    arrange(standtree, mm) %>% 
  # below line is to correct rings starting above 1. If they start at 1, has no    effect
    mutate(ring = ring - (min(ring) - 1))
  
  df.1 <- df %>% 
    group_by(standtree, ring) %>% 
    summarise(start = min(mm), end = max(mm)) %>% 
    mutate(ring = seq(from = 1, to = (max(ring)-min(ring)) + 1))
  
  df2 <- df %>% 
    left_join(df.1, by = c("standtree", "ring")) %>% 
    group_by(standtree, ring) %>% 
    mutate(cutoff = min(density) + (max(density)-min(density)) * .8,
           latewood = ifelse(density > cutoff, 1, 0))
  
  # Step 2 and 3: Identify and assign new ring latewood boundaries and ring        boundaries, based on first occurrence of latewood
  
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

# This step is why this needs to be done on a scan-by-scan basis. If I re-wrote this to work with lists, I could probably correct that.
  df3 <- df2 %>% 
    ungroup() %>% 
    group_by(standtree) %>% 
    mutate(newring = as.numeric(cut(df2$mm, breaks = c(0,as.vector(
      df.2$boundary2))))) 
  
  df3 <- df3 %>% 
    left_join(df.2, by = c("newring" = "ring")) %>% 
    select(-standtree.y) %>% 
    rename(standtree = standtree.x)
  
  # Step 4: Join boundaries onto original dataframe by "newring". 
  
  df4 <- df3 %>% 
    ungroup() %>% 
    mutate(newEorL = ifelse(mm < boundary1, "EW", "LW")) %>% 
    select(standtree, mm, newring, density, cutoff, newEorL) %>% 
    as.data.frame()
  
  corScan <<- df4
  }