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
calc_chem_tox_hits <- function(chem_tox_dat, chem_conc_dat, sites_detect, sites_hit) {
  
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

merge_top_chems_hits <- function(tox_hits, wq_hits, top_chems) {
  dat <- left_join(top_chems, wq_hits) %>%
    left_join(tox_hits)
  
  return(dat)
}

summarize_chems <- function(file_name, chem_vals, chem_crosswalk, chem_info) {
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

calc_parent_tox_hits <- function(parent_sums){
  hits_ear <- parent_sums %>%
    filter(type == 'p_sumval' & measure_type %in% 'ear') %>%
    filter(sumval > 0.001) %>%
    group_by(parent_pesticide) %>%
    summarize(n_hits = n(), 
              n_hits_sites = length(unique(site)),
              n_hits_months = length(unique(lubridate::month(date))),
              median_sumval = median(sumval),
              max_sumval = max(sumval)) %>%
    mutate(measure_type = 'ear')
  
  hits_bench <- parent_sums %>%
    filter(type == 'p_sumval' & measure_type %in% 'bench') %>%
    filter(sumval > 0.01) %>%
    group_by(parent_pesticide) %>%
    summarize(n_hits = n(), 
              n_hits_sites = length(unique(site)),
              n_hits_months = length(unique(lubridate::month(date))),
              median_sumval = median(sumval),
              max_sumval = max(sumval)) %>%
    mutate(measure_type = 'bench')
}
