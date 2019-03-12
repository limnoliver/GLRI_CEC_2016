# summary fig of data availability
library(dplyr)
library(tidyr)
library(remake)
dat <- make('merged_dat', remake_file = '20_merge_data.yml')
owc <- make('owc')
head(dat)

sites <- make('sites') %>%
  select(site_no, shortName) %>%
  unique()

dat_all <- owc %>%
  mutate(source = 'owc_s4433') %>%
  bind_rows(dat) %>%
  left_join(sites, by = c('SiteID' = 'site_no'))

  
dat_summary <- dat_all %>%
  group_by(SiteID, shortName, source) %>%
  summarize(n_dates = length(unique(sample_dt))) %>%
  spread(key = source, value = n_dates)

dat_month_summary <- dat_all %>%
  mutate(month = lubridate::month(sample_dt, label = TRUE, abbr = TRUE)) %>%
  group_by(SiteID, shortName, source, month) %>%
  summarize(n_samples = length(unique(sample_dt))) %>%
  spread(key = source, value = n_samples)

dat_month_summary[is.na(dat_month_summary)] <- 0

glyphosate_month <- dat_month_summary %>%
  select(-neonic, -pesticides_s2437, -owc_s4433) %>%
  spread(key = month, value = glyphosate)

neonic_month <- dat_month_summary %>%
  select(-glyphosate, -pesticides_s2437, -owc_s4433) %>%
  spread(key = month, value = neonic)

pesticides_month <- dat_month_summary %>%
  select(-glyphosate, -neonic, -owc_s4433) %>%
  spread(key = month, value = pesticides_s2437)

owc_month <- dat_month_summary %>%
  select(-glyphosate, -neonic, -pesticides_s2437) %>%
  spread(key = month, value = owc_s4433)

write.csv(dat_summary, '30_reports/GLRI_2016_pesticides_owc_summary.csv', row.names = F)
write.csv(glyphosate_month, '30_reports/GLRI_2016_glyphosate_monthly_summary.csv', row.names = F)
write.csv(neonic_month, '30_reports/GLRI_2016_neonic_monthly_summary.csv', row.names = F)
write.csv(pesticides_month, '30_reports/GLRI_2016_pesticides_monthly_summary.csv', row.names = F)
write.csv(owc_month, '30_reports/GLRI_2016_owc_monthly_summary.csv', row.names = F)

head(dat_summary)

# QA data

rep_dat <- filter(dat_all, samp_type_cd == 7)

rep_summary <- rep_dat %>%
  group_by(SiteID, shortName, source) %>%
  summarize(n_dates = length(unique(sample_dt))) %>%
  spread(key = source, value = n_dates)
