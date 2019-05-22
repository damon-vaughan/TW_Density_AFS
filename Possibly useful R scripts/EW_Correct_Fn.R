# Correct full profiles so that EW density equals min EW density until the mm position where min EW is reached. Basically fixes curve rings and early compression wood.

EW.fix <- function(x){
  xlist <- split(x, f = x$newring)
  
  x.minEW <- x %>% 
    group_by(newring) %>% 
    summarise(minmm = min(mm),
              minew = which.min(density) * .01 - .01,
              minEW = minmm + minew) %>% 
    select(newring, minEW)
  
  xlist2 <- lapply(xlist, function(x){left_join(x, x.minEW, by = "newring")})
  
  x2 <- lapply(xlist2, function(x){df <- x %>% 
    mutate(density = ifelse(mm < minEW, min(density), density))}) %>% 
    rbindlist()
}

