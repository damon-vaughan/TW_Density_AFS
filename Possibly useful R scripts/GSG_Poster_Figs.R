Fig3a <- ggplot(filter(climWide, grow.yr >= 1940 & grow.yr <= 2016)) +
  geom_line(aes(x = grow.yr, y = PDSIavg)) +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("PDSI") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3b <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = RD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("RD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3c <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = EWD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2, alpha = .5) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             alpha = .5, color = "red") +
  ylab("EWD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(1940, 2010, by = 10))

Fig3aGrob <- ggplotGrob(Fig3a)
Fig3bGrob <- ggplotGrob(Fig3b)
Fig3cGrob <- ggplotGrob(Fig3c)

Fig3Grob <- rbind(Fig3aGrob, Fig3bGrob, Fig3cGrob, size = "first")

Fig3 <- grid.arrange(Fig3Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/PosterFig2.tiff",
       device = "tiff", plot = Fig3,
       width = 28, height = 14, units = "cm")
