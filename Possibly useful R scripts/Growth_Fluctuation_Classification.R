df3 <- df2 %>% 
  filter(standtree == "TW-5-366-2") 

ggplot(df3) + 
  geom_smooth(aes(x = year, y = RW), method = "gam", formula = y ~ s(x, k = 5)) +
  geom_smooth(aes(x = year, y = RW), method = "lm")
              
              +
  geom_hline(yintercept = mean(df3$RW))

library(mgcv)
mod1 <- mgcv::gam(RW ~ s(year, k = 5), data = df3)
mod2 <- lm(RW ~ year, data = df3)

newdat <- data.frame(year = df3$year,
                     gam = predict(mod1), 
                     lm = predict(mod2))

newdat <- newdat %>% 
  mutate(change = ifelse(gam < lm, -1, 1),
         lag = lag(change, 1),
         flag = change + lag) 
row.names(newdat) <- df3$year

xvec <- rownames(newdat)[which(newdat$flag == 0)]

ggplot(newdat) +
  geom_bar(aes(x = year, y = change), stat = "identity")
