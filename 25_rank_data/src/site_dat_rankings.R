# create site-specific info
get_all_samples <- function(chem_data) {
  unique_samples <- select(chem_data, site = SiteID, date = `Sample Date`) %>%
    distinct()
  
  return(unique_samples)
}
get_nchems_sites <- function(chemicalSummary_conc, all_samples) {
  nchems <- chemicalSummary_conc %>%
    group_by(site, date) %>%
    summarize(n_detected = n())
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    filter()
  
  nchems <- full_join(nchems, samples) %>%
    mutate(n_detected = ifelse(is.na(n_detected), 0, n_detected)) %>%
    arrange(site, date)
  
  # how many chemicals are unique to this sample within this site?
  unique_samples <- select(chemicalSummary_conc, site, date) %>%
    distinct()
  
  unique_samples$n_new <- NA
  for (i in 1:nrow(unique_samples)) {
    temp_site <- unique_samples$site[i]
    temp_date <- as.Date(unique_samples$date[i])
    
    other_site_dat <- filter(chemicalSummary_conc, site %in% temp_site) %>%
      filter(date != temp_date)
    other_chems <- unique(other_site_dat$CAS)
    
    sample_dat <- filter(chemicalSummary_conc, site %in% temp_site) %>%
      filter(date == temp_date)
    sample_chems <- unique(sample_dat$CAS)
    
    unique_samples$n_new[i] <- length(which(!(sample_chems %in% other_chems)))
  }
  
  nchems <- left_join(nchems, unique_samples) %>%
    mutate(n_new = ifelse(is.na(n_new), 0, n_new))
    
  
  nchem_ranks <- nchems %>%
    group_by(site) %>%
    summarize(med_ndetect = median(n_detected), 
              max_ndetect = max(n_detected),
              med_ndetect_new = median(n_new), 
              max_ndetect_new = max(n_new)) %>%
    mutate(med_ndetect_rank = rank(-med_ndetect),
           max_ndetect_rank = rank(-max_ndetect), 
           med_ndetect_new_rank = rank(-med_ndetect_new),
           max_ndetect_new_rank = rank(-max_ndetect_new)) %>%
    select(site, contains('rank'))
  
  nchems_wide <- left_join(nchems, nchem_ranks)
  
  nchem_ndetect <- nchems_wide %>%
    select(-contains('new')) %>%
    rename(value = n_detected,
           median_rank = med_ndetect_rank, 
           max_rank = max_ndetect_rank) %>%
    mutate(variable = 'n_detected')
  
  nchem_ndetect_new <- nchems_wide %>%
    select(site, date, contains('new')) %>%
    rename(value = n_new,
           median_rank = med_ndetect_new_rank, 
           max_rank = max_ndetect_new_rank) %>%
    mutate(variable = 'n_detected_new')
  
  nchems <- bind_rows(nchem_ndetect, nchem_ndetect_new)
  
  return(nchems)
}

get_conc_sites <- function(chemicalSummary_conc, all_samples) {
  site_conc <- chemicalSummary_conc %>%
    group_by(site, date) %>%
    summarize(sum_conc = sum(EAR), 
              mean_conc = mean(EAR)) %>%
    ungroup()
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    filter()
  
  site_conc <- full_join(site_conc, samples) %>%
    mutate(sum_conc = ifelse(is.na(sum_conc), 0, sum_conc),
           mean_conc = ifelse(is.na(mean_conc), 0, mean_conc))
  
  site_conc_rank <- group_by(site_conc, site) %>%
    summarize(med_sum_conc = median(sum_conc),
           max_sum_conc = max(sum_conc),
           med_mean_conc = median(mean_conc),
           max_mean_conc = max(mean_conc)) %>%
    mutate(med_sum_conc_rank = rank(-med_sum_conc),
           max_sum_conc_rank = rank(-max_sum_conc),
           med_mean_conc_rank = rank(-med_mean_conc),
           max_mean_conc_rank = rank(-max_mean_conc)) %>%
    select(site, contains('rank'))
  
  
  site_conc_wide <- left_join(site_conc, site_conc_rank)
  
  sum_conc <- select(site_conc_wide, -contains('mean_conc')) %>%
    rename(value = sum_conc, 
           median_rank = med_sum_conc_rank,
           max_rank = max_sum_conc_rank) %>%
    mutate(variable = 'sum_conc')
  
  mean_conc <- select(site_conc_wide, -contains('sum_conc')) %>%
    rename(value = mean_conc, 
           median_rank = med_mean_conc_rank,
           max_rank = max_mean_conc_rank) %>%
    mutate(variable = 'mean_conc')
  
  site_conc <- bind_rows(sum_conc, mean_conc)
                           
  return(site_conc)
  
}

get_ear_sites <- function(chemicalSummary, all_samples) {
  site_ear <- chemicalSummary %>%
    group_by(site, date, CAS) %>%
    summarize(sum_ear = sum(EAR)) %>%
    ungroup() %>%
    group_by(site, date) %>%
    summarize(mean_sumEAR = mean(sum_ear),
              max_sumEAR = max(sum_ear))
  
  # this likely underestimates effects, but summing likely overestimates
  # mean might not make sense -- having a "low" value doesn't somehow offset a 
  # high value. 
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    filter()
  
  site_ear <- full_join(site_ear, samples) %>%
    mutate(mean_sumEAR = ifelse(is.na(mean_sumEAR), 0, mean_sumEAR),
           max_sumEAR = ifelse(is.na(max_sumEAR), 0, max_sumEAR)) %>%
    ungroup()
  
  site_ear_rank <- group_by(site_ear, site) %>%
    summarize(med_mean_sumEAR = median(mean_sumEAR),
              max_mean_sumEAR = max(mean_sumEAR),
              med_max_sumEAR = median(max_sumEAR),
              max_max_sumEAR = max(max_sumEAR)) %>%
    mutate(med_mean_sumEAR_rank = rank(-med_mean_sumEAR),
           max_mean_sumEAR_rank = rank(-max_mean_sumEAR),
           med_max_sumEAR_rank = rank(-med_max_sumEAR),
           max_max_sumEAR_rank = rank(-max_max_sumEAR)) %>%
    select(site, contains('rank'))
  
  site_ear_wide <- left_join(site_ear, site_ear_rank)
  
  site_mean_ear <- site_ear_wide %>%
    select(-contains('max_sumEAR')) %>%
    rename(value = mean_sumEAR,
           median_rank = med_mean_sumEAR_rank,
           max_rank = max_mean_sumEAR_rank) %>%
    mutate(variable = 'mean_sumEAR')
  
  site_max_ear <- site_ear_wide %>%
    select(-contains('mean_sumEAR')) %>%
    rename(value = max_sumEAR,
           median_rank = med_max_sumEAR_rank,
           max_rank = max_max_sumEAR_rank) %>%
    mutate(variable = 'max_sumEAR')
  
  site_ear <- bind_rows(site_mean_ear, site_max_ear)
  
  return(site_ear)
  
}

get_bench_sites <- function(chemicalSummary_bench, all_samples) {
  site_bench <- chemicalSummary_bench %>%
    group_by(site, date) %>%
    summarize(mean_bench = mean(EAR), 
              max_bench = max(EAR))
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    filter()
  
  site_bench <- full_join(site_bench, samples) %>%
    mutate(mean_bench = ifelse(is.na(mean_bench), 0, mean_bench),
           max_bench = ifelse(is.na(max_bench), 0, max_bench)) %>%
    ungroup()
  
  site_bench_rank <- group_by(site_bench, site) %>% 
    summarize(med_mean_bench = median(mean_bench),
           max_mean_bench = max(mean_bench),
           med_max_bench= median(max_bench),
           max_max_bench = max(max_bench)) %>%
    mutate(med_mean_bench_rank = rank(-med_mean_bench),
           max_mean_bench_rank = rank(-max_mean_bench),
           med_max_bench_rank = rank(-med_max_bench),
           max_max_bench_rank = rank(-max_max_bench)) %>%
    select(site, contains('rank'))
  
  site_bench_long <- left_join(site_bench, site_bench_rank)
  
  site_bench_mean <- site_bench_long %>%
    select(-contains('max_bench')) %>%
    rename(value = mean_bench,
           median_rank = med_mean_bench_rank,
           max_rank = max_mean_bench_rank) %>%
    mutate(variable = 'mean_bench')
  
  site_bench_max <- site_bench_long %>%
    select(-contains('mean_bench')) %>%
    rename(value = max_bench,
           median_rank = med_max_bench_rank,
           max_rank = max_max_bench_rank) %>%
    mutate(variable = 'max_bench')
  
  site_bench <- bind_rows(site_bench_mean, site_bench_max)
  
  return(site_bench)
  
}

merge_site_rankings <- function(nchems, conc, ear, bench) {
  rankings <- left_join(nchems, conc) %>%
    left_join(ear) %>%
    left_join(bench)
  
  vals <- rankings %>%
    select(-contains('rank')) %>%
    gather(key = 'rank_')
  
  return(rankings)
}

calc_avg_rankings <- function(site_rankings, sites) {
  unique_rank <- select(site_rankings, site, median_rank, max_rank, variable) %>%
    distinct()
  
  median_rank <- select(unique_rank, -max_rank) %>%
    spread(key = variable, value = median_rank)
  
  max_rank <- select(unique_rank, -median_rank) %>%
    spread(key = variable, value = max_rank)
  
  library(matrixStats)
  
  ranks <- left_join(median_rank, max_rank, by = 'site') %>%
    ungroup() %>%
    mutate(avg_rank = rowMeans(select(.,`max_bench.x`:`sum_conc.y`)),
           sum_rank = rowSums(select(.,`max_bench.x`:`sum_conc.y`)),
            sd_rank = rowSds(as.matrix(select(.,`max_bench.x`:`sum_conc.y`))))
  
  ranks_site <- select(ranks, site, avg_rank, sd_rank) %>%
    left_join(sites, by = c('site' = 'site_no')) %>%
    mutate(disturbance_index = `Ag..crops` + Urban, 
           natural_index = Forest + `Water..wetland`) %>%
    mutate(dist_diff_index = disturbance_index - natural_index)
  
  return(ranks_site)
  
}

plot_avg_rankings <- function(outfile, ranking_dat, sites) {
  
  dat <- left_join(ranking_dat, sites, by = c('site' = 'SiteID'))
 
   p <- ggplot(dat, aes(y = avg_rank, x = disturbance_index, label = `Short Name`)) +
    geom_point(aes(color = `Dominant.land.use.`), size = 2) + 
    geom_errorbar(aes(ymin = avg_rank - sd_rank, 
                      ymax = avg_rank + sd_rank, 
                      color = `Dominant.land.use.`), width = 1.5) +
     #geom_label(check_overlap = T, hjust = -0.1, nudge_x = 0.3) +
     ggrepel::geom_label_repel(data = dat, aes(color = `Dominant.land.use.`), show.legend = F, size = 2) +
    theme_bw() +
    labs(x = "Watershed Disturbance Index (% Crop + % Urban)", 
         y = 'Average (sd) rank across metrics\n(1 = most impacted by pesticides)',
         color = "Dominant land use") +
    scale_x_continuous(limits = c(0,100))
  
  ggsave(outfile, p, height = 4, width = 8)
}
