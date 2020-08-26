library(ggplot2)
library(dplyr)


chem_bench <- make('chemicalSummary_bench_deg_meto') %>%
  group_by(site, date, CAS, chnm, parent_pesticide) %>%
  summarize(maxBench = max(EAR)) %>%
  group_by(site, date, parent_pesticide) %>%
  summarize(summaxBench = sum(maxBench))

test <- ungroup(chem_bench) %>%
  left_join(graph_dat)
#levels(chem_bench$chnm)[86] <- '2,4-D'

parents <- make('chem_master') %>%
  filter(!grepl('Deg', Class)) %>%
  select(parent_pesticide, Class) %>%
  distinct()

drop_dates <- make('maumee_exclude')
graph_dat <- make('parent_sums') %>%
  left_join(parents) %>%
  mutate(day = lubridate::yday(date)) %>%
  filter(type == 'p_d_sumval' & measure_type == 'bench') %>%
  filter(sumval >= 0.1) %>%
  group_by(site, day, date, Class) %>%
  summarize(n_hits = n(),
            max_bench = max(sumval)) %>%
  ungroup() %>%
  filter(!(site %in% '04193500' & as.Date(date) %in% as.Date(drop_dates)))

sites <- make('sites') %>%
  select(site = site_no, dominant_lu = Dominant.land.use.)

graph_dat <- left_join(graph_dat, sites) %>%
  add_row(day = c(1, 365), Class = rep('Fungicide', 2), dominant_lu = rep('Crops', 2)) %>%
  filter(!dominant_lu == 'Wetland')

p <- ggplot(data = graph_dat, aes(x = day, y = max_bench)) +
  geom_hline(yintercept = c(0.1, 1, 10), color = 'gray50') +
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
  labs(x = '', y = 'Max TQchem for each site-date-class', color = 'Dominant Land Use', size = '# Chems w/hits')

ggsave('figures/ms_figures/month_class_hits_bench.png', p, height = 6, width = 8)
