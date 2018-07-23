# occurance calculations by chemical

chem_data <- make('chem_data')
chem_info <- make('chem_info')
sum_conc <- make('sum_conc')
dat <- group_by(chem_data, pCode, CAS) %>%
  summarise(count_detected = length(Value[!(remark_cd %in% "<")]), 
            count_all = n(),
            count_sites = length(unique(SiteID[!(remark_cd %in% "<")])),
            mean = mean(Value[!(remark_cd %in% "<")]), 
            median = median(Value[!(remark_cd %in% "<")]),
            min = min(Value),
            max = max(Value), 
            sd = sd(Value[!(remark_cd %in% "<")])) %>%
  mutate(detection_freq = count_detected/count_all)

dat <- left_join(dat, chem_info, by = 'CAS')

dat_insect <- filter(dat, Class %in% c('Deg - Insecticide', 'Insecticide'))
dat_fung <- filter(dat, Class %in% c('Fungicide', 'Deg - Fungicide'))

# summary stats from site perspective
site_dat <- filter(chem_data, !(remark_cd %in% "<")) %>% 
  group_by(SiteID) %>%
  summarise(unique_chem = length(unique(pCode)))

site_dat_avg <- sum_conc %>%
  group_by(SiteID) %>%
  summarize(n_chem_mean = mean(n_detected),
            sum_conc_mean = mean(sum_conc_detect))
                     
