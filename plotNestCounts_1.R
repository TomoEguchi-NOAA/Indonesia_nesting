#plotNestCounts


rm(list=ls())
source('Dc_Indonesia_nesting_fcns.R')

saveFig <- F

data.0 <- read.csv('data/NestCounts_22Jan2018.csv')

# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0 %>% mutate(begin.date = as.Date(paste(Year_begin,
                                             Month_begin,
                                             '01', sep = "-"),
                                       format = "%Y-%m-%d")) %>%
  mutate(end.day = if_else(Month_end == 1 | Month_end == 3 | Month_end == 5 |
                             Month_end == 7 | Month_end == 8 | Month_end == 10 |
                             Month_end == 12, 31,
                           if_else(Month_end == 2,
                                   if_else(leap_year(Year_end) == TRUE,
                                           29, 28), 30))) %>%
  mutate(end.date = as.Date(paste(Year_end, Month_end, end.day,
                                  sep = "-"),
                            format = "%Y-%m-%d")) %>%
  select(begin.date, end.date, J_M, Com_JM, W, Com_W) %>%
  select(begin.date, end.date, J_M, W, Com_JM, Com_W) %>%
  mutate(cumu_days = as.numeric(begin.date - begin.date[1])) %>%
  gather(., location, count, J_M:W) -> data.1

# fix one end date:
# data.1[grep("UntilOct18", data.1$Com_W)[2],
#        'end.date'] <- '2016-10-18'
data.1$num.days = as.numeric(data.1$end.date - data.1$begin.date + 1)

data.1.JM <- filter(data.1, location == 'J_M') %>%
  select(., begin.date, end.date, count, num.days) %>%
  mutate(Year_begin = year(begin.date)) %>%
  group_by(., Year_begin) %>%
  summarise(., sum.nests = sum(count), sum.time = sum(num.days)) %>%
  mutate(., nests.per.time = sum.nests/sum.time) %>%
  na.omit()

fit.JM <- lm(log(nests.per.time) ~ Year_begin, data = data.1.JM)
predict.JM <- predict(fit.JM, se.fit = T)
data.1.JM$predicted <- predict.JM$fit
data.1.JM$predict.SE <- predict.JM$se.fit

data.1.W <- filter(data.1, location == 'W') %>%
  select(., begin.date, end.date, count, num.days) %>%
  mutate(Year_begin = year(begin.date)) %>%
  group_by(., Year_begin) %>%
  summarise(., sum.nests = sum(count), sum.time = sum(num.days)) %>%
  mutate(., nests.per.time = sum.nests/sum.time) %>%
  na.omit()

fit.W <- lm(log(nests.per.time) ~ Year_begin, data = data.1.W)
predict.W <- predict(fit.W, se.fit = T)
data.1.W$predicted <- predict.W$fit
data.1.W$predict.SE <- predict.W$se.fit

plot.1 <- ggplot() +
  geom_point(data = data.1.JM,
             aes(x = Year_begin, y = nests.per.time),
                 color = 'blue', size = 3) +
  geom_line(data = data.1.JM,
            aes(x = Year_begin,
                y = exp(predicted)),
            color = 'blue', size = 2) +
  geom_line(data = data.1.JM,
            aes(x = Year_begin,
                y = exp(predicted + 2 * predict.SE)),
            color = 'blue', size = 1.3) +
  geom_line(data = data.1.JM,
            aes(x = Year_begin,
                y = exp(predicted - 2 * predict.SE)),
            color = 'blue', size = 1.3) +

  geom_point(data = data.1.W,
             aes(x = Year_begin, y = nests.per.time),
             color = 'red', size = 3) +
  geom_line(data = data.1.W,
            aes(x = Year_begin, y = exp(predicted)),
            color = 'red', size = 2) +
  geom_line(data = data.1.W,
            aes(x = Year_begin,
                y = exp(predicted + 2 * predict.SE)),
            color = 'red', size = 1.3) +
  geom_line(data = data.1.W,
            aes(x = Year_begin,
                y = exp(predicted - 2 * predict.SE)),
            color = 'red', size = 1.3) +

  labs(x = '', y = '# nests/sample duration (days)')  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)


# ggsave(plot.1,
#        filename = 'figures/Indonesia_trend.png',
#        dpi = 600)

# data.2.JM <- filter(data.1.JM, Year_begin < 2012)
# fit.JM.2 <- lm(log(nests.per.time) ~ Year_begin, data = data.2.JM)
#
# data.2.W <- filter(data.1.W, Year_begin < 2012)
# fit.W.2 <- lm(log(nests.per.time) ~ Year_begin, data = data.2.W)

data.1.2005 <- filter(data.1, begin.date >= as.Date("2005-01-01",
                                                    format = "%Y-%m-%d")) %>%
  mutate(cumu_days = as.numeric(begin.date - as.Date("2005-01-01",
                                          format = "%Y-%m-%d"))) %>%
  select(., -starts_with('Com'))
  #na.omit()


plot.2 <- ggplot(data = data.1.2005) +
  geom_point(aes(x = begin.date,
                 y = count/num.days,
                 color = location),
             size = 2) +
  geom_line(aes(x = begin.date,
                 y = count/num.days,
                 color = location),
             size = 1) +

  labs(x = '', y = '# nests/day')  +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Jamursba-Medi", "Wermon"),
                     name = "Location") +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90,
                                   hjust = 0.5,
                                   vjust = 0.5),
        legend.position = c(0.7, 0.9),
        legend.title.align = 0.5)


if (saveFig)
  ggsave(plot.2,
         filename = 'figures/Indonesia_trend_since2005.png',
         dpi = 600)
