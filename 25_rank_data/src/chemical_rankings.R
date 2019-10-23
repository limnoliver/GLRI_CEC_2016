# get chemicals that are in >50% of months at >50% of sites
library(lubridate)
select_top_chems <- function(chem_dat) {
  count_by_site_month <- mutate(chem_dat, month = month(date)) %>%
    group_by(site, CAS, chnm) %>%
    summarize(month_count = length(unique(month)))
  
  chems_count_sites <- group_by(count_by_site_month, CAS, chnm) %>%
    summarize(n_sites_2months = length(which(month_count >= 2)),
              n_sites = n())
  
  top_chems <- filter(chems_count_sites, n_sites_2months >= 8)
  return(top_chems)
  
  top_chems_attr <- left_join(top_chems, distinct(select(chem_dat, CAS, Class)))
}
calc_chem_tox_hits <- function(chem_tox_dat, chem_conc_dat) {
  
  site_hits_tox <- chem_tox_dat %>%
    group_by(site, date, CAS, chnm) %>%
    summarize(sumEAR = sum(EAR)) %>%
    filter(sumEAR > 0.001) %>%
    ungroup() %>%
    group_by(CAS, chnm) %>%
    summarize(n_hits = n(),
              n_site_hits = length(unique(site))) %>%
    filter(n_site_hits >= sites_hit)
  
  # site_samples <- reduced_dat %>%
  #   group_by(SiteID, sample_dt) %>%
  #   summarize(n_chems = n()) %>%
  #   group_by(SiteID) %>%
  #   summarize(n_samples = n()) %>%
  #   rename(site = SiteID)
  # 
  # hits_per_site <- chem_tox_dat %>%
  #   group_by(site, date, CAS, chnm) %>%
  #   summarize(sumEAR = sum(EAR)) %>%
  #   filter(sumEAR > 0.001) %>%
  #   ungroup() %>%
  #   group_by(site, CAS, chnm) %>%
  #   summarize(n_hits_per_site = n()) %>%
  #   left_join(site_samples) %>%
  #   mutate(prop_sample_hits = n_hits_per_site/n_samples) %>%
  #   group_by(CAS, chnm) %>%
  #   summarize(mean_prop_sample_hits = mean(prop_sample_hits))
    
  
  site_detects <- chem_conc_dat %>%
    group_by(CAS, chnm) %>%
    summarize(n_sites = length(unique(site))) %>%
    filter(n_sites >= sites_detect)
  
  top_chems <- inner_join(site_hits_tox, site_detects, by = c('CAS', 'chnm'))
  return(top_chems)
}

calc_chem_tox_metrics <- function(chem_dat) {
 
  hits_tox <- chem_dat %>%
    group_by(site, date, CAS, chnm) %>%
    summarize(sumEAR = sum(EAR)) %>%
    group_by(CAS, chnm) %>%
    summarize(tox_n_detect = n(), 
              n_hits_tox = length(which(sumEAR > 0.01)), 
              n_sites_tox_hits = length(unique(site[which(sumEAR > 0.01)])),
              mean_sumEAR = mean(sumEAR),
              max_sumEAR = max(sumEAR))
  
  return(hits_tox)
}

calc_chem_wq_metrics <- function(chem_dat) {

  hits_wq <- chem_dat %>%
    group_by(site, date, CAS, chnm) %>%
    summarize(max_bench = max(EAR)) %>%
    group_by(CAS, chnm) %>%
    summarize(wq_n_detect = n(), 
              n_hits_wq = length(which(max_bench > 0.1)), 
              n_sites_hits = length(unique(site[which(max_bench > 0.1)])),
              mean_bench = mean(max_bench),
              max_bench = max(max_bench))
  
  return(hits_wq)
  
}

cal_chem_mix_metrics <- function(mix_dat) {
  
}

merge_top_chems_hits <- function(tox_hits, wq_hits, top_chems) {
  dat <- left_join(top_chems, wq_hits) %>%
    left_join(tox_hits)
  
  return(dat)
}

summarize_chems <- function(file_name, chem_vals, chem_crosswalk, chem_info, conc_dat) {
  sample_count <- all_chems %>%
    group_by(CAS, pCode) %>%
    summarize(n_samples = n(),
              n_sites = length(unique(SiteID)))
  
  n_detect <- all_chems %>%
    filter(!remark_cd %in% '<') %>%
    group_by(CAS) %>%
    summarize(n_detect = n(),
              n_detect_sites = length(unique(SiteID)))
  
  detects <- left_join(sample_count, n_detect) %>%
    mutate(n_detect = ifelse(is.na(n_detect), 0, n_detect),
           n_detect_sites = ifelse(is.na(n_detect_sites), 0, n_detect_sites),
           perc_detect = round(100*(n_detect/n_samples), 0),
           perc_detect_site = round(100*(n_detect_sites/n_sites), 0)) %>%
    mutate(detection_column = paste0(perc_detect, ' (', n_detect, '/', n_samples, ')'),
           detection_site_column = paste0(perc_detect_site, ' (', n_detect_sites, '/', n_sites, ')')) %>%
    select(CAS, pCode, detection_column, detection_site_column) %>%
    left_join(select(chem_crosswalk, CAS, compound, parent_pesticide)) %>%
    left_join(chem_info)
  
  write.csv(detects, file_name, row.names = FALSE)
    
}

summarize_chem_meta <- function(file_name, chem_vals, chem_crosswalk, chem_info_all, chems_missing_toxcast, chems_missing_bench, neonic) {
  
  sample_count <- chem_vals %>%
    group_by(CAS, pCode) %>%
    summarize(n_samples = n(),
              n_sites = length(unique(SiteID)))
  
  n_detect <- chem_vals %>%
    filter(!remark_cd %in% '<') %>%
    group_by(CAS) %>%
    summarize(n_detect = n(),
              n_detect_sites = length(unique(SiteID))) %>%
    mutate(n_detect = ifelse(is.na(n_detect), 0, n_detect),
           n_detect_sites = ifelse(is.na(n_detect_sites), 0, n_detect_sites))

  meta <- left_join(sample_count, chem_info_all, by = 'CAS') %>%
    left_join(select(chem_crosswalk, CAS, compound, parent_pesticide)) %>%
    mutate(`In toxCast` = ifelse(CAS %in% chems_missing_toxcast, '', 'x'),
           `In benchmarks` = ifelse(CAS %in% chems_missing_bench, '', 'x')) %>%
    left_join(n_detect, by = 'CAS')
  
  meta$parent_pesticide[is.na(meta$parent_pesticide)] <- meta$`Chemical Name`[is.na(meta$parent_pesticide)]
  
  meta <- meta %>%
    arrange(parent_pesticide) %>%
    select(`Chemical Name`, CAS, `USGS parameter code` = pCode, 
           Class, `Parent compound` = parent_pesticide, `In toxCast`, 
           `In benchmarks`, `N samples measured` = n_samples, `N samples detected` = n_detect, 
           `N sites measured` = n_sites, `N sites detected` = n_detect_sites)
  
  
  write.csv(meta, file_name, row.names = FALSE)
  
}

calc_parent_tox_hits <- function(parent_sums, mixtures){
  hits_ear <- parent_sums %>%
    filter(type == 'p_d_sumval' & measure_type %in% 'ear') %>%
    #filter(sumval > 0.001) %>%
    group_by(parent_pesticide) %>%
    summarize(ear_hit_prob = round(length(which(sumval > 0.001))/n(),2), 
              ear_hit_sites_prob = round(length(unique(site[sumval > 0.001]))/length(unique(site)), 2),
              ear_hit_months = length(unique(lubridate::month(date[sumval > 0.001]))),
              ear_median_sumval = round(median(sumval), 4),
              ear_max_sumval = round(max(sumval), 4))
  
  hits_bench <- parent_sums %>%
    filter(type == 'p_d_sumval' & measure_type %in% 'bench') %>%
    #filter(sumval > 0.01) %>%
    group_by(parent_pesticide) %>%
    summarize(bench_hit_prob = round(length(which(sumval > 0.01))/n(), 2), 
              bench_hit_sites_prob = round(length(unique(site[sumval > 0.01]))/length(unique(site)), 2),
              bench_hit_months = length(unique(lubridate::month(date[sumval > 0.01]))),
              bench_median_sumval = round(median(sumval), 3),
              bench_max_sumval = round(max(sumval), 3))
  
  hits_mixes <- mixtures %>%
    rename(mixes_hit_n = times_in_mixes, mixes_hit_sites_n = n_sites, 
           mixes_hit_endpoints_n = n_endpoints, mixes_median_contribution = contribution_median) %>%
    mutate(mixes_median_contribution_prob = round(mixes_median_contribution/100, 2)) %>%
    select(-mixes_median_contribution)
  
  detections <- parent_sums %>%
    filter(type == 'p_d_sumval' & measure_type %in% 'conc') %>%
    group_by(parent_pesticide) %>%
    summarize(detect_prob = round(length(which(sumval > 0))/n(), 2),
              detect_sites_prob = round(length(unique(site[sumval > 0]))/length(unique(site)),2),
              detect_months = length(unique(lubridate::month(date[sumval > 0]))))
  
  parent_metrics <- left_join(detections, hits_mixes, by = 'parent_pesticide') %>%
    left_join(hits_bench, by = 'parent_pesticide') %>%
    left_join(hits_ear, by = 'parent_pesticide')
  
  return(parent_metrics)
    
}

filter_parents <- function(parent_metrics, f_detect_prob, f_detect_sites_prob, f_detect_months, f_hit_prob, f_hit_sites_prob) {
  top_parents <- parent_metrics %>%
    filter(detect_prob >= f_detect_prob |
             detect_sites_prob >= f_detect_sites_prob |
             detect_months >= f_detect_months |
             bench_hit_prob >= f_hit_prob |
             ear_hit_prob >= f_hit_prob |
             ear_hit_sites_prob >= f_hit_sites_prob |
             ear_hit_sites_prob >= f_hit_sites_prob
             )
}
