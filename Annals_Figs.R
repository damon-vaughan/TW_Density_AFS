
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Axis titles at 12 pt, axis text at 10 pt. Same with legends

# Fig 1 -------------------------------------------------------------------

Fig1a <- ggplot(PRCP) +
  geom_line(aes(x = grow.yr, y = PRCPtot)) +
  geom_smooth(aes(x = grow.yr, y = PRCPtot), method = "gam", se = FALSE, 
              formula = y ~ s(x, k = 5)) +
  ylab("Total Precipitation (mm)") +
  xlab("Water Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))  +
  scale_x_continuous(breaks = seq(1920, 2020, by = 20))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig1a.eps",
       plot = Fig1a, device = "eps", 
       width = 129, height = 80, units = "mm")

Fig1b <- ggplot(Normals) +
  geom_line(aes(x = Month, y = PRCPmonth)) +
  ylab("Total Precipitation (mm)") +
  xlab("Month") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(1,12, by = 1))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig1b.eps",
       plot = Fig1b, device = "eps", 
       width = 43, height = 40, units = "mm")

# Fig 2 -------------------------------------------------------------------

Fig2 <- ggplot(TW_Growth) +
  stat_summary(aes(x = year, y = BA, color = factor(GSL)), 
               fun.y = mean, geom = "line", size = 0.5) +
  geom_hline(yintercept = c(6.9, 13.8, 18.4, 23, 27.6, 34.4), 
             color = cbPalette[c(2,3,4,6,7,8)],
             linetype = "dashed", size = 0.5) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  ylab(expression(paste("Basal area ", "(m" ^ "2 ","ha" ^ "-1", ")"))) +
  labs(color = bquote(atop(GSL~phantom(),(m^2~ha^-1)))) +
  scale_x_continuous(breaks = seq(1960, 2020, by = 20)) +
  scale_color_manual(values = cbPalette[c(2,3,4,6,7,8)]) +
  guides(color = guide_legend(reverse = TRUE))

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig2.eps",
       device = "eps", plot = Fig2,
       width = 84, height = 55, units = "mm")


# Fig 3 -------------------------------------------------------------------
# Had to take out alpha adjustment to save as eps

Fig3a <- ggplot(filter(climWide, grow.yr >= 1940 & grow.yr <= 2016)) +
  geom_line(aes(x = grow.yr, y = PDSIavg)) +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             color = "red") +
  ylab("PDSI") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3b <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = RD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             color = "red") +
  ylab("RD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3c <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = EWD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             color = "red") +
  ylab("EWD") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3d <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = LWP), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2) +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             color = "red") +
  ylab("LWP") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

Fig3e <- ggplot(filter(df3, year >= 1940 & year <= 2016)) +
  stat_summary(aes(x = year, y = LWD), 
               fun.y = function(x) mean(x), na.rm = T, geom = "line") +
  geom_vline(xintercept = c(1951, 1956, 1963, 1977, 1989, 1996, 2002), 
             linetype = 2) +
  ylab("LWD") +
  geom_vline(xintercept = c(1962, 1972, 1982, 1992, 2003),
             color = "red") +
  scale_x_continuous(breaks = seq(1940, 2010, by = 10)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank())

Fig3aGrob <- ggplotGrob(Fig3a)
Fig3bGrob <- ggplotGrob(Fig3b)
Fig3cGrob <- ggplotGrob(Fig3c)
Fig3dGrob <- ggplotGrob(Fig3d)
Fig3eGrob <- ggplotGrob(Fig3e)

Fig3Grob <- rbind(Fig3aGrob, Fig3bGrob, Fig3cGrob, 
                  Fig3dGrob, Fig3eGrob, size = "first")

Fig3 <- grid.arrange(Fig3Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig3.eps",
       device = "eps", plot = Fig3,
       width = 174, height = 174, units = "mm")


# Fig 4 -------------------------------------------------------------------
# TAke off axis.line = black, because i'm leaving the panel borders on

Fig4a <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9), 
                aes(x = factor(GSL), y = BAImean, color = factor(Sample.HT))) +
  geom_boxplot(size = .75, notch = F, width = .6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  ylab(expression(paste("BAI ", "(cm" ^ "2", ")"))) +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

Fig4b <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9), 
                aes(x = factor(GSL), y = RDmean, color = factor(Sample.HT))) +
  geom_boxplot(size = .75, notch = F, width = .6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  ylab(expression(paste("RD ", "(kg m" ^ "-3", ")"))) +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

Fig4c <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9),  
                aes(x = factor(GSL), y = LWPmean, color = factor(Sample.HT))) +
  geom_boxplot(size = .75, notch = F, width = .6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 15)),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  xlab(expression(paste("GSL ", "(m" ^ "2 ","ha" ^ "-1", ")"))) +
  ylab("LWP") +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

Fig4aGrob <- ggplotGrob(Fig4a)
Fig4bGrob <- ggplotGrob(Fig4b)
Fig4cGrob <- ggplotGrob(Fig4c)
Fig4Grob <- rbind(Fig4aGrob, Fig4bGrob, Fig4cGrob, size = "first")
Fig4 <- grid.arrange(Fig4Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig4a.eps", plot = Fig4,
       device = "eps", 
       width = 174, height = 174, units = "mm")

Fig4legend <- ggplot(filter(dfsum, Sample.HT == 1.4 | Sample.HT == 2.4 | Sample.HT == 4.9),   
                     aes(x = factor(GSL), y = BAImean)) +
  geom_boxplot(aes(color = factor(Sample.HT)), size = .75) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  labs(color = "Sample\nheight (m)") +
  scale_colour_manual(values = colorschemes$SteppedSequential.5[c(6,13,18)])

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig4Legend.eps", 
       plot = Fig4legend,
       device = "eps", 
       width = 174, height = 174, units = "mm")

col2rgb(colorschemes$SteppedSequential.5[c(6,13,18)])

# Fig 5 -------------------------------------------------------------------

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Fig5a <- ggplot(filter(df62, Sample.HT == 1.4), 
                aes(x = year, y = BAI, color = factor(GSL))) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), se = F, size = 0.5) +
  ylab(expression(paste("BAI ", "(cm" ^ "2", ")"))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 12, margin = margin(r = 15))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(1963, 2016), breaks = seq(1970, 2010, by = 10)) +
  # scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  scale_color_manual(values = cbPalette[c(2,3,4,6,7,8)]) +
  theme(plot.margin = unit(c(t = .5, r = .5, b = .5, l = .5), "cm"))

Fig5b <- ggplot(filter(df62, Sample.HT == 1.4), 
                aes(x = year, y = RD, color = factor(GSL))) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), se = F, size = 0.5) +
  ylab(expression(paste("RD ", "(kg m" ^ "-3", ")"))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 12, margin = margin(r = 15))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(1963, 2016), breaks = seq(1970, 2010, by = 10)) +
  # scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  scale_color_manual(values = cbPalette[c(2,3,4,6,7,8)]) +
  theme(plot.margin = unit(c(t = .5, r = .5, b = .5, l = .5), "cm"))

Fig5c <- ggplot(filter(df62, Sample.HT == 1.4),  
                aes(x = year, y = LWP, color = factor(GSL))) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), se = F, size = 0.5) +
  ylab("LWP") +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 12, margin = margin(r = 15))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(1963, 2016), breaks = seq(1970, 2010, by = 20)) +
  # scale_color_manual(values = colorschemes$Categorical.12[c(1,6,8,10,11,12)]) +
  scale_color_manual(values = cbPalette[c(2,3,4,6,7,8)]) +
  theme(plot.margin = unit(c(t = .5, r = .5, b = .5, l = .5), "cm"))

Fig5aGrob <- ggplotGrob(Fig5a)
Fig5bGrob <- ggplotGrob(Fig5b)
Fig5cGrob <- ggplotGrob(Fig5c)
Fig5Grob <- rbind(Fig5aGrob, Fig5bGrob, Fig5cGrob, size = "first")
Fig5 <- grid.arrange(Fig5Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig5.eps", plot = Fig5,
       device = "eps", 
       width = 75, height = 110, units = "mm")

Fig5Legend <- ggplot(df62) +
  geom_line(aes(x = year, y = BAI, color = factor(GSL)), size = .5) +
  theme_bw() +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  labs(color = bquote(atop(GSL~phantom(),(m^2~ha^-1)))) +
  scale_color_manual(values = cbPalette[c(2,3,4,6,7,8)])

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig5_Legend.eps", plot = Fig5Legend,
       device = "eps", 
       width = 75, height = 110, units = "mm")

col2rgb(colorschemes$Categorical.12[c(1,6,8,10,11,12)])


# Fig 6 -------------------------------------------------------------------

RQ2.PRCP <- read.csv("RQ2_Results_PRCP.csv", header = T)
RQ2.TAVG <- read.csv("RQ2_Results_TAVG.csv", header = T)

RQ2.PRCP$group <- factor(RQ2.PRCP$group, levels = c("Low", "Mid", "High"))
RQ2.TAVG$group <- factor(RQ2.TAVG$group, levels = c("Low", "Mid", "High"))

# Top row plots

Fig6a <- ggplot(filter(RQ2.PRCP, Response == "RD"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  ylab("PRCP\nresponse coef.") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  ggtitle("RD") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 12)) +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig6b <- ggplot(filter(RQ2.PRCP, Response == "LWP"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none") +
  ggtitle("LWP") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 12)) +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig6c <- ggplot(filter(RQ2.PRCP, Response == "EWD"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none") +
  ggtitle("EWD") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 12)) +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig6d <- ggplot(filter(RQ2.PRCP, Response == "MXD"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.PRCP$quarter[1:5]) +
  scale_y_continuous(limits = c(-.5,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none") +
  ggtitle("MXD") + 
  theme(plot.title = element_text(margin = margin(b = -5), 
                                  hjust = .5, size = 12)) +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig6aGrob <- ggplotGrob(Fig6a)
Fig6bGrob <- ggplotGrob(Fig6b)
Fig6cGrob <- ggplotGrob(Fig6c)
Fig6dGrob <- ggplotGrob(Fig6d)

Fig6_1Grob <- cbind(Fig6aGrob, Fig6bGrob, Fig6cGrob, Fig6dGrob, size = "first")

Fig6.1 <- grid.arrange(Fig6_1Grob)

# Bottom row plots

Fig7a <- ggplot(filter(RQ2.TAVG, Response == "RD"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  ylab("TAVG\nresponse coef.") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig7b <- ggplot(filter(RQ2.TAVG, Response == "LWP"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig7c <- ggplot(filter(RQ2.TAVG, Response == "EWD"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig7d <- ggplot(filter(RQ2.TAVG, Response == "MXD"),
                aes(x = quarter, y = coef, ymin = ci_lower, 
                    ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(position = position_dodge(width = .5), width = 0.4, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(limits = RQ2.TAVG$quarter[1:5]) +
  scale_y_continuous(limits = c(-.41,.5)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cbPalette[c(3,4,8)])

Fig7aGrob <- ggplotGrob(Fig7a)
Fig7bGrob <- ggplotGrob(Fig7b)
Fig7cGrob <- ggplotGrob(Fig7c)
Fig7dGrob <- ggplotGrob(Fig7d)

Fig6_2Grob <- cbind(Fig7aGrob, Fig7bGrob, Fig7cGrob, Fig7dGrob, size = "first")

Fig6Grob <- rbind(Fig6_1Grob, Fig6_2Grob, size = "first")

Fig6 <- grid.arrange(Fig6Grob)

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig6.eps", plot = Fig6,
       device = "eps", 
       width = 160, height = 160, units = "mm")

# Legend material

Fig6_2Legend <- ggplot(RQ2.TAVG,
                       aes(x = quarter, y = coef, ymin = ci_lower, 
                           ymax = ci_upper, color = factor(group))) +
  geom_point(position = position_dodge(width = .75)) +
  geom_errorbar(position = position_dodge(width = .75), width = 0.4, size = .5) +
  ylab("Temperature\nresponse coef.") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  labs(color = "GSL group") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = cbPalette[c(3,4,8)])

ggsave("C:/Users/vaug8/OneDrive/WQA/Documentation/TW_Paper_2019/Images/Fig6_2Legend.eps", plot = Fig6_2Legend,
       device = "eps", 
       width = 174, height = 174, units = "mm")

col2rgb(cbPalette[c(3,4,8)])
