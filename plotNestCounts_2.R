#plotNestCounts


rm(list=ls())
source('Dc_Indonesia_nesting_fcns.R')

saveFig <- T

data.raw <- read.csv('data/adjusted number of nests.csv')

data.0 <- gather(data.raw, Location, Count, Jamursba_Medi, Wermon)

plot.1 <- ggplot(data = data.0) +
  geom_point(aes(x = Year, y = Count,
                 color = Location), size = 3) +
  geom_line(aes(x = Year, y = Count,
                 color = Location), size = 1) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Jamursba-Medi", "Wermon"),
                     name = "Location") +
  scale_x_continuous(breaks = seq(1980, 2020, by=5)) +
  scale_y_continuous(breaks = seq(0, 16000, by=2000),
                     limits = c(0, 15000)) +
  labs(x = '', y = '# nests')  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.8, 0.8),
        legend.title.align = 0.5)

plot.1

plot.1.log <- ggplot(data = data.0) +
  geom_point(aes(x = Year, y = log(Count),
                 color = Location), size = 3) +
  geom_line(aes(x = Year, y = log(Count),
                color = Location), size = 1) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Jamursba-Medi", "Wermon"),
                     name = "Location") +
  scale_x_continuous(breaks = seq(1980, 2020, by=5)) +
  #scale_y_continuous(breaks = seq(0, 16000, by=2000),
  #                   limits = c(0, 15000)) +
  labs(x = '', y = 'log(# nests)')  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.8, 0.8),
        legend.title.align = 0.5)

plot.1.log

data.1 <- filter(data.0, Year > 2000)
plot.2 <- ggplot(data = data.1) +
  geom_point(aes(x = Year, y = Count,
                 color = Location), size = 3) +
  geom_line(aes(x = Year, y = Count,
                color = Location), size = 1) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Jamursba-Medi", "Wermon"),
                     name = "Location") +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
  scale_y_continuous(breaks = c(0, 2000, 4000),
                     limits = c(0, 5000)) +
  labs(x = '', y = '# nests')  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.8, 0.8),
        legend.title.align = 0.5)

plot.2

plot.2.log <- ggplot(data = data.1) +
  geom_point(aes(x = Year, y = log(Count),
                 color = Location), size = 3) +
  geom_line(aes(x = Year, y = log(Count),
                color = Location), size = 1) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Jamursba-Medi", "Wermon"),
                     name = "Location") +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
  #scale_y_continuous(breaks = c(0, 2000, 4000),
  #                   limits = c(0, 5000)) +
  labs(x = '', y = 'log(# nests)')  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.8, 0.8),
        legend.title.align = 0.5)

plot.2.log

# plot.2 <- ggplot(data = data.1.2005) +
#   geom_point(aes(x = begin.date,
#                  y = count/num.days,
#                  color = location),
#              size = 2) +
#   geom_line(aes(x = begin.date,
#                  y = count/num.days,
#                  color = location),
#              size = 1) +
#
#   labs(x = '', y = '# nests/day')  +
#   scale_color_manual(values = c("blue", "red"),
#                      labels = c("Jamursba-Medi", "Wermon"),
#                      name = "Location") +
#   scale_x_date(date_breaks = "1 year",
#                date_labels = "%Y") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         axis.text.x = element_text(angle = 90,
#                                    hjust = 0.5,
#                                    vjust = 0.5),
#         legend.position = c(0.7, 0.9),
#         legend.title.align = 0.5)

if (saveFig){
  ggsave(plot.1,
         filename = 'figures/Indonesia_nest_trend.png',
         dpi = 600, height = 6, width = 8, units = "in")
  ggsave(plot.2,
         filename = 'figures/Indonesia_nest_trend_since2001.png',
         dpi = 600, height = 6, width = 8, units = "in")

  ggsave(plot.1.log,
         filename = "figures/Indoensia_nest_trend_log.png",
         dpi = 600, height = 6, width = 8, units = "in")

  ggsave(plot.2.log,
         filename = "figures/Indoensia_nest_trend_log_since2001.png",
         dpi = 600, height = 6, width = 8, units = "in")

}
