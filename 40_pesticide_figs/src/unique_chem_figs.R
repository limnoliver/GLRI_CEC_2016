get_n_unique_chems <- function(reduced_dat, sites) {
  
  n_unique_time <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("SiteID", "sample_dt", "unique_on_date"))
  site_ids <- unique(reduced_dat$SiteID)
  site_num_names <- select(sites, site_no, shortName, `Dominant.land.use.`)
  
  for (i in 1:length(site_ids)) {
    temp_dat <- filter(reduced_dat, SiteID == site_ids[i]) %>%
      filter(!(remark_cd %in% '<')) %>%
      select(SiteID, sample_dt, pCode) %>%
      arrange(sample_dt)
    
    dates <- unique(temp_dat$sample_dt)
    
    
    for (j in 1:length(dates)) {
      
      # calculate number of new chems
      ordered_dat <- filter(temp_dat, sample_dt <= dates[j]) %>%
        summarize(unique_on_date = length(unique(pCode))) %>%
        mutate(SiteID = site_ids[i], sample_dt = dates[j])
      
      n_unique_time <- bind_rows(n_unique_time, ordered_dat)
      
    }
  }
  
  n_unique_time <- group_by(n_unique_time, SiteID) %>%
    mutate(n_added_on_date = c(NA, diff(unique_on_date))) %>%
    mutate(n_added_on_date = ifelse(is.na(n_added_on_date), unique_on_date, n_added_on_date)) %>%
    mutate(sample_no = 1:n()) %>%
    left_join(site_num_names, by = c('SiteID' = 'site_no'))
  
  return(n_unique_time)
  
}


plot_unique_chems <- function(target_name, unique_chems, type){
  
  if (type == 'cumulative') {
    p <- ggplot(n_unique_time, aes(x = sample_dt, y = unique_on_date)) +
      geom_line(aes(group = shortName, color = `Dominant.land.use.`)) +
      labs(x = '', y = 'Unique chemicals detected', color = "Dominant land use") +
      theme_bw()
  }
  
  if (type == 'new') {
    # n_after_first <- filter(n_unique_time, sample_no >1)
    # ggplot(n_after_first, aes(sample_dt, y = n_added_on_date)) +
    #   geom_line(aes(group = shortName, color = shortName)) +
    #   facet_wrap(~shortName, ncol = 2, scales = 'free_y') +
    #   theme_bw()
    
    p <- ggplot(n_unique_time, aes(sample_dt, y = n_added_on_date)) +
      geom_line(aes(group = shortName, color = `Dominant.land.use.`)) +
      facet_wrap(~shortName, ncol = 2, scales = 'free_y') +
      theme_bw() +
      labs(y = 'N unique chemicals added', x = '', color = 'Dominant land use')
  }
  
  ggsave(filename = target_name, p, height = ifelse(type == 'new', 12, 6), width = ifelse(type == 'new', 8, 6))
  
  
  
}



