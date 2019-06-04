test <- data.frame(chem = sample(c('a', 'b', 'c', 'd', 'e', 'f'), 100, replace = TRUE),
                   value = sample(seq(0.001, 1, by = 0.001), 100, replace = TRUE),
                   month = sample(c(1:12), 100, replace = TRUE))

test <- mutate(test, class = case_when(
  chem %in% c('a', 'b') ~ 'class1',
  chem %in% c('c', 'd') ~ 'class2',
  chem %in% c('e', 'f') ~ 'class3'
)) %>%
  add_row(month = 13)
  
head(test)
library(ggplot2)
library(dplyr)
ggplot(data = test, aes(x = month, y = value)) +
  geom_point(aes(color = class, shape = chem), alpha = 0.5) +
  #geom_jitter() +
  coord_polar(theta = 'x', start = 0) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,4,7,10), labels = c('Jan', 'Apr', 'Jul', 'Oct'))

chem_ear <- make('chemicalSummary_deg_meto')
#levels(chem_ear$chnm)[86] <- '2,4-D'

parents <- chem_ear %>%
  filter(!grepl('deg', Class, ignore.case = TRUE)) %>%
  select(Class, parent_pesticide) %>%
  distinct()

chem_sums <- chem_ear %>%
  group_by(site, date, parent_pesticide) %>%
  summarize(sumEAR = sum(EAR)) %>%
  left_join(parents)

graph_dat <- make('parent_sums') %>%
  left_join(make('parent_class')) %>%
  mutate(day = lubridate::yday(date)) %>%
  filter(type == 'p_sumval' & measure_type == 'ear') %>%
  filter(sumval >= 0.001) %>%
  group_by(site, day, Class) %>%
  summarize(n_hits = n(),
            max_ear = max(sumval)) %>%
  ungroup()

sites <- make('sites') %>%
  select(site = site_no, dominant_lu = Dominant.land.use.)

graph_dat <- left_join(graph_dat, sites) %>%
  add_row(day = c(1, 365), Class = rep('Herbicide', 2), dominant_lu = rep('Crops', 2)) %>%
  filter(!dominant_lu == 'Wetland')

p <- ggplot(data = graph_dat, aes(x = day, y = max_ear)) +
  geom_hline(yintercept = c(0.001, 0.01, 0.1, 1), color = 'gray50') +
  geom_vline(xintercept = c(1,32,61, 92,122,153,183,214,245,275, 306, 336), color = 'gray80') +
  geom_point(aes(color = site, size = n_hits), alpha = 0.5, shape = 16) +
  facet_grid(rows = vars(dominant_lu), cols = vars(Class)) +
  #geom_jitter() +
  coord_polar(theta = 'x', start = 0) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,32,61, 92,122,153,183,214,245,275, 306, 336), labels = c('Jan','','Mar', '', 'May', '', 'Jul', '', 'Sep', '', 'Nov', '')) +
  scale_y_log10() +
  theme(panel.grid  = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold')) +
  labs(x = '', y = 'Max EAR for each site-date-class', color = 'Dominant Land Use', size = '# Chems w/hits')

ggsave('40_pesticide_figs/month_class_hits.png', p, height = 6, width = 8)
