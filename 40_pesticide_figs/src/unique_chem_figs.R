# number of unique chemicals through time by site
reduced_dat <- make('reduced_dat')
sites <- unique(reduced_dat$SiteID)

n_unique_time <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("SiteID", "sample_dt", "unique_on_date"))

for (i in 1:length(sites)) {
  temp_dat <- filter(reduced_dat, SiteID == sites[i]) %>%
    filter(!(remark_cd %in% '<')) %>%
    select(SiteID, sample_dt, pCode) %>%
    arrange(sample_dt)
  
  dates <- unique(temp_dat$sample_dt)

  
  for (j in 1:length(dates)) {
    
    # calculate number of new chems
    ordered_dat <- filter(temp_dat, sample_dt <= dates[j]) %>%
      summarize(unique_on_date = length(unique(pCode))) %>%
      mutate(SiteID = sites[i], sample_dt = dates[j])
    
    n_unique_time <- bind_rows(n_unique_time, ordered_dat)
      
  }
}

n_unique_time <- group_by(n_unique_time, SiteID) %>%
  mutate(n_added_on_date = c(NA, diff(unique_on_date))) %>%
  mutate(n_added_on_date = ifelse(is.na(n_added_on_date), unique_on_date, n_added_on_date)) %>%
  mutate(sample_no = 1:n())

ggplot(n_unique_time, aes(x = sample_dt, y = unique_on_date)) +
  geom_line(aes(group = SiteID, color = SiteID)) +
  labs(x = '', y = 'Unique chemicals detected') +
  theme_bw()

n_after_first <- filter(n_unique_time, sample_no >1)
ggplot(n_after_first, aes(sample_dt, y = n_added_on_date)) +
  geom_line(aes(group = SiteID, color = SiteID)) +
  facet_wrap(~SiteID, ncol = 2, scales = 'free_y') +
  theme_bw()
