# create site-specific info
get_all_samples <- function(chem_data) {
  unique_samples <- select(chem_data, site = SiteID, date = `Sample Date`) %>%
    distinct()
  
  return(unique_samples)
}
get_nchems_sites <- function(chemicalSummary_conc, all_samples) {
  nchems <- chemicalSummary_conc %>%
    filter(EAR > 0) %>%
    group_by(site, date) %>%
    summarize(n_detected = n())
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date))
  
  nchems <- full_join(nchems, samples) %>%
    mutate(n_detected = ifelse(is.na(n_detected), 0, n_detected)) %>%
    arrange(site, date)
    
  nchem_unique_by_site <- chemicalSummary_conc %>%
    filter(EAR > 0) %>%
    group_by(site) %>%
    summarize(n_detect_all = length(unique(chnm)))
    
  
  nchem_metrics <- nchems %>%
    group_by(site) %>%
    summarize(median_ndetect = median(n_detected), 
              max_ndetect = max(n_detected)) %>%
    left_join(nchem_unique_by_site)
  
  nchem_metrics_long <- nchem_metrics %>%
    tidyr::gather(key = 'metric', value = 'metric_value', -site)
  
  nchem_ranks_long <- nchem_metrics_long %>%
    group_by(metric) %>%
    mutate(rank_value = rank(-metric_value),
           relative_value = metric_value/max(metric_value),
           metric_type = 'Occurance')
  

  # nchem_ndetect <- nchems_wide %>%
  #   select(-contains('new')) %>%
  #   rename(value = n_detected,
  #          median_rank = med_ndetect_rank, 
  #          max_rank = max_ndetect_rank) %>%
  #   mutate(variable = 'n_detected')
  # 
  # nchem_ndetect_new <- nchems_wide %>%
  #   select(site, date, contains('new')) %>%
  #   rename(value = n_new,
  #          median_rank = med_ndetect_new_rank, 
  #          max_rank = max_ndetect_new_rank) %>%
  #   mutate(variable = 'n_detected_new')
  
  #nchems <- bind_rows(nchem_ndetect, nchem_ndetect_new)
  
  return(nchem_ranks_long)
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
  site_sum_ear <- chemicalSummary %>%
    #group_by(site, date, CAS) %>%
    group_by(site, date, parent_pesticide) %>%
    summarize(sum_ear = sum(EAR))
  
  site_ear_summary <- site_sum_ear %>%
    group_by(site, date) %>%
    summarize(max_sumEAR = max(sum_ear),
              n_hits = length(which(sum_ear > 0.001))) %>%
    right_join(mutate(all_samples, date = as.POSIXct(date))) %>%
    mutate(max_sumEAR = ifelse(is.na(max_sumEAR), 0, max_sumEAR),
           n_hits = ifelse(is.na(n_hits), 0, n_hits)) %>%
    group_by(site) %>%
    summarize(max_max_sumEAR = round(max(max_sumEAR),4),
              med_max_sumEAR = round(median(max_sumEAR),4),
              med_hits_per_sample = median(n_hits),
              n_month_hits = length(unique(lubridate::month(date)[n_hits > 0]))) %>%
    select(site, max_max_sumEAR, med_max_sumEAR, med_hits_per_sample, n_month_hits) 
  
  site_ear_metric_long <- site_ear_summary %>%
    tidyr::gather(key = 'metric', value = 'metric_value', -site)
  
  site_ear_rank_long <- site_ear_metric_long %>%
    group_by(metric) %>%
    mutate(rank_value = rank(-metric_value),
           relative_value = metric_value/max(metric_value)) %>%
    mutate(metric_type = 'EAR')
  
  return(site_ear_rank_long)
  
}

get_bench_sites <- function(chemicalSummary_bench, all_samples) {
  
  # first, find lowest toxicity unit/highest TQ for each chemical
  # (conservative approach)
  chem_max_bench <- chemicalSummary_bench %>%
    group_by(site, date, CAS, parent_pesticide) %>%
    summarize(max_bench = max(EAR))
  
  # now, sum the TQs for each parent compound
  parent_sum_bench <- chem_max_bench %>%
    group_by(site, date, parent_pesticide) %>%
    summarize(sum_max_bench = sum(max_bench))
  
  site_sum_ear <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    left_join(parent_sum_bench) %>%
    mutate(sum_max_bench = ifelse(is.na(sum_max_bench), 0, sum_max_bench))
  
  sites_bench_vals <- group_by(parent_sum_bench, site, date) %>%
    summarize(max_sumBench = max(sum_max_bench),
              n_hits = length(which(sum_max_bench > 0.1))) %>%
    right_join(mutate(all_samples, date = as.POSIXct(date))) %>%
    mutate(max_sumBench = ifelse(is.na(max_sumBench), 0, max_sumBench),
           n_hits = ifelse(is.na(n_hits), 0, n_hits)) %>%
    group_by(site) %>%
    summarize(max_max_sumBench = round(max(max_sumBench),4),
              med_max_sumBench = round(median(max_sumBench),4),
              med_hits_per_sample = median(n_hits),
              n_month_hits = length(unique(lubridate::month(date)[n_hits > 0]))) %>%
    select(site, max_max_sumBench, med_max_sumBench, med_hits_per_sample, n_month_hits) 
  
  
  site_bench_metric_long <- sites_bench_vals %>%
    tidyr::gather(key = 'metric', value = 'metric_value', -site)
  
  site_bench_rank_long <- site_bench_metric_long %>%
    group_by(metric) %>%
    mutate(rank_value = rank(-metric_value),
           relative_value = metric_value/max(metric_value)) %>%
    mutate(metric_type = 'Bench')
  
  
  return(site_bench_rank_long)
  
}

get_mix_sites <- function(mixes) {
  mixes_long <- mixes %>%
    select(-shortName) %>%
    rename(Hits = mix_hits_per_sample,
           `max EARmix` = max_EARmix, 
           `Months w/hits` = mix_hit_n_months,
           `Endpoints w/hits` = mix_hit_n_endpoints) %>%
    tidyr::gather(key = 'metric', value = 'metric_value', -site) %>%
    group_by(metric) %>%
    mutate(rank_value = rank(-metric_value),
           relative_value = metric_value/max(metric_value)) %>%
    mutate(metric_type = 'EARmix')
  
  return(mixes_long)
    
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
