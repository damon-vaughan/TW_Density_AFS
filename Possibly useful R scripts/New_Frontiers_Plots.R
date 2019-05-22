ggplot(df2) +
  # geom_point(aes(x=year, y=ringDens), col = "gray50") +
  geom_smooth(aes(x=year, y=ringDens, color = factor(Sample.pos)),
              method = "loess",
              # se = FALSE,
              size = 1) +
  # geom_vline(xintercept = 1962) +
  ylab(expression(paste("Ring density  ", "(kg ", "m"^"-3",")"))) +
  xlab("Cambial Age")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 24)) +
  scale_color_discrete("Height (ft)")

ggsave("Images/NewFrontiers/MeanRingDens.tiff", 
       device = "tiff", 
       width = 22, height = 14, units = "cm")

#############################
# Plots of the components- EW dens, LW dens, and widths

# Ring density
p1 <- ggplot(sumProfile.Int) +
  geom_smooth(aes(x = year, y = ringDens), method = "gam",
              formula = y ~ s(x, k = 5)) +
  # stat_summary(aes(x = year, y = ringDens, color = factor(HT)), fun.y = mean,
  # geom = "line") +
  ylab(expression(paste("Ring density ", "(kg ", "m"^"-3",")"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(1930, 2010, by = 20))

# LW proportion
p2 <- ggplot(df2) +
  geom_smooth(aes(x = year, y = LW.prop), method = "gam",
              formula = y ~ s(x, k = 5)) +
  stat_summary(aes(x = year, y = LW.prop, color = factor(HT)), fun.y = mean,
               geom = "line") +
  ylab("LW proportion") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(1930, 2010, by = 20))

# EW density
p3 <- ggplot(df2) +
  geom_smooth(aes(x = year, y = EW.dens), method = "gam",
              formula = y ~ s(x, k = 5)) +
  stat_summary(aes(x = year, y = EW.dens, color = factor(HT)), fun.y = mean,
               geom = "line") +
  ylab(expression(paste("EW density ", "(kg ", "m"^"-3",")"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(1930, 2010, by = 20))

# LW density
p4 <- ggplot(df2) +
  geom_smooth(aes(x = year, y = LW.dens), method = "gam",
              formula = y ~ s(x, k = 5)) +
  stat_summary(aes(x = year, y = LW.dens, color = factor(HT)), fun.y = mean,
               geom = "line") +
  ylab(expression(paste("LW density ", "(kg ", "m"^"-3",")"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(1930, 2010, by = 20))

library(STA578)
multiplot(p1,p3,p2,p4, ncol = 2)
# Export at 1000/600

#############
# EW/LW Width and Ring Width
ggplot(df2) +
  geom_smooth(aes(x = year, y = LW.width), color = "red",
              method = "gam", formula = y ~ s(x, k = 8)) +
  annotate(geom="text", x=2000, y=.4, label="Latewood Width",
           color="red", size = 6) +
  geom_smooth(aes(x = year, y = EW.width), color = "blue", 
              method = "gam", formula = y ~ s(x, k = 8)) +
  annotate(geom="text", x=1980, y=1.7, label="Earlywood width",
           color="blue", size = 6) +
  geom_smooth(aes(x = year, y = RW), color = "black",
              method = "gam", formula = y ~ s(x, k = 8)) +
  annotate(geom="text", x=2005, y=2, label="Total ring width",
           color="black", size = 6) +
  ylab("Width (mm)") +
  xlab("Year") +
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16))

ggsave("Images/NewFrontiers/RW_components.tiff", 
       device = "tiff", 
       width = 22, height = 14, units = "cm")

#################################
# Last minute: diff heights over year

ggplot(df2) +
  stat_summary(aes(x=year, y = ringDens, color = factor(HT)), fun.y = mean,
               geom = "line") +
  ylab(expression(paste("Ring density  ", "(kg ", "m"^"-3",")"))) +
  xlab("Year") +
  theme(axis.title.y = element_text(size = 20,
                                    margin = margin(
                                      t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 24)) +
  scale_color_discrete("Height (ft)") +
  # geom_vline(xintercept = c(2002, 1996, 1989, 1977, 1974,
  # 1971, 1963, 1956, 1951, 1947)) +
  # xlim(1920,2017) +
  scale_x_continuous(breaks = seq(1930, 2010, by = 20))

ggsave("Images/NewFrontiers/DiffHT_vsYear.tiff", 
       device = "tiff", 
       width = 22, height = 14, units = "cm")

##########################
# Multiplot of yearly ring density and weather fluctuations
p1 <- ggplot(w.sum) +
  geom_line(aes(x = grow.yr, y = rain), color = "blue") +
  theme(axis.title.y = element_text(size = 20,
                                    margin = margin(
                                      t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Total rainfall (in)") +
  xlab("Year") +
  geom_vline(xintercept = c(2002, 1996, 1989, 1977, 1974,
                            1971, 1963, 1956, 1951, 1947)) +
  xlim(1920,2017) 

p2 <- ggplot(df2) +
  stat_summary(aes(x=year, y = ringDens), fun.y = mean,
               geom = "line", color = "Red") +
  ylab(expression(paste("Ring density  ", "(kg ", "m"^"-3",")"))) +
  xlab("Year") +
  theme(axis.title.y = element_text(size = 20,
                                    margin = margin(
                                      t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16)) +
  geom_vline(xintercept = c(2002, 1996, 1989, 1977, 1974,
                            1971, 1963, 1956, 1951, 1947)) +
  # xlim(1920,2017) +
  scale_x_continuous(breaks = seq(1930, 2010, by = 20))

ggsave(plot = p2, "Images/NewFrontiers/DensYear.tiff", 
       device = "tiff", 
       width = 30, height = 14, units = "cm")

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size = "first")
grid.newpage()
grid.draw(g)
