library(tidyverse)

dat <- fullProfile 

dat2 <- dat %>% 
  group_by(standtree, newring) %>% 
  mutate(EorL.EF = ifelse(density >= cutoff, "L", "E"),
         lagEorL.EF = lag(EorL.EF),
         EtoLmark = ifelse(EorL.EF == "L" & lagEorL.EF == "E", 1, 0)) %>% 
  na.omit()

dat3 <- dat2 %>% 
  group_by(standtree, newring) %>% 
  mutate(position = cumsum(EtoLmark),
         error1 = ifelse(density >= cutoff, abs(density-cutoff), 0),
         error2 = ifelse(density < cutoff, abs(density-cutoff), 0)) 

dat4 <- dat3 %>% 
  group_by(standtree, newring, position) %>% 
  summarise(Emagnitude = sum(error1),
            Elength = length(error1),
            Lmagnitude = sum(error2),
            Llength = length(error2)) %>% 
  filter(position > 0) %>% 
  mutate(EF1 = Elength^3 * Emagnitude,
         EF2 = Llength^3 * Lmagnitude) %>% 
  ungroup() %>% 
  select(standtree, newring, position, EF1, EF2)

dat_nest <- dat4 %>% 
  nest(-standtree, -newring)

EF.func <- function(x) {
  
  newdf <- data.frame(matrix(ncol = 2)) %>% 
    rename("position" = X1, "EF" = X2)
  
  max.p <- max(x$position)
  for (m in 1:max.p) {
    newdf[m, 1] = m
  }
  
  for (i in 2:max.p) {
    newdf[1, 2] <- sum(x[2:max.p, 2])
    j = i - 1
    k = i + 1
    newdf[i, 2] <- sum(x[1:j, 3], x[k:max.p, 2])
    newdf[max.p, 2] <- sum(x[1:(max.p - 1), 3])
  }
  return(which.min(newdf$EF))
  # return(newdf)
}

min.EF <- dat_nest %>%
  mutate(min.EF = map(data, EF.func)) %>%
  select(-data) %>% 
  unnest(min.EF)

str(min.EF)

dat5 <- dat3 %>% 
  left_join(min.EF, by = c("standtree", "newring")) %>% 
  ungroup()

dat6 <- dat5 %>% 
  mutate(newLmark = ifelse(position == min.EF & EtoLmark == 1, 1, 0))

dat7 <- dat6 %>% 
  group_by(standtree, newring) %>% 
  mutate(newEorL.EF = cumsum(newLmark))

dat8 <- dat7 %>% 
  mutate(newEorL.EF2 = ifelse(newEorL.EF == 1, "L", "E"))


