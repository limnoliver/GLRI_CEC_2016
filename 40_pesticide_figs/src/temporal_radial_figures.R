library(ggplot2)
library(dplyr)


chem_ear <- make('chemicalSummary_deg_meto')
#levels(chem_ear$chnm)[86] <- '2,4-D'

parents <- chem_ear %>%
  filter(!grepl('deg', Class, ignore.case = TRUE)) %>%
  select(Class, parent_pesticide) %>%
  distinct()

drop_dates <- make('maumee_exclude')


chem_sums <- chem_ear %>%
  group_by(site, date, parent_pesticide) %>%
  summarize(sumEAR = sum(EAR)) %>%
  left_join(parents)

graph_dat <- make('parent_sums') %>%
  left_join(make('parent_class')) %>%
  mutate(day = lubridate::yday(date)) %>%
  filter(type == 'p_d_sumval' & measure_type == 'ear') %>%
  filter(sumval >= 0.001) %>%
  filter(!(site %in% '04193500' & as.Date(date) %in% as.Date(drop_dates))) %>%
  group_by(site, day, Class) %>%
  summarize(n_hits = n(),
            max_ear = max(sumval)) %>%
  ungroup()

sites <- make('sites') %>%
  select(site = site_no, dominant_lu = Dominant.land.use.)

sites2 <- make('sites') %>%
  select(site = site_no, shortName, dominant_lu = Dominant.land.use.)
graph_dat <- left_join(graph_dat, sites) %>%
  add_row(day = c(1, 365), Class = rep('Herbicide', 2), dominant_lu = rep('Crops', 2)) %>%
  filter(!dominant_lu == 'Wetland')

p <- ggplot(data = graph_dat, aes(x = day, y = max_ear)) +
  geom_hline(yintercept = c(0.001, 0.01, 0.1, 1), color = 'gray50') +
  geom_vline(xintercept = c(1,32,61, 92,122,153,183,214,245,275, 306, 336), color = 'gray80') +
  geom_point(aes(size = n_hits), alpha = 0.5, shape = 16, color = 'darkred') +
  facet_grid(rows = vars(dominant_lu), cols = vars(Class)) +
  #geom_jitter() +
  coord_polar(theta = 'x', start = 0) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,32,61, 92,122,153,183,214,245,275, 306, 336), labels = c('Jan','','Mar', '', 'May', '', 'Jul', '', 'Sep', '', 'Nov', '')) +
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid  = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold')) +
  labs(x = '', y = 'Max EAR for each site-date-class', color = 'Dominant Land Use', size = '# Chems w/hits')

ggsave('figures/ms_figures/month_class_hits.png', p, height = 6, width = 8)

head(graph_dat)
ggplot(graph_dat, aes(x = day, y = n_hits)) +
  geom_point(aes(group = site, color = site)) +
  geom_line(aes(group = site, color = site)) +
  facet_wrap(~dominant_lu)

graph_dat <- make('parent_sums') %>%
  left_join(make('parent_class')) %>%
  filter(type == 'p_d_sumval' & measure_type == 'ear') %>%
  filter(sumval >= 0.001) %>%
  filter(!(site %in% '04193500' & as.Date(date) %in% as.Date(drop_dates))) %>%
  mutate(month = lubridate::month(date)) %>%
  group_by(site) %>%
  summarize(n_month_hits = length(unique(month))) %>%
  left_join(sites2)

            