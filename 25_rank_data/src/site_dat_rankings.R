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
    
    other_site_dat <- filter(chemicalSummary, site %in% temp_site) %>%
      filter(date != temp_date)
    other_chems <- unique(other_site_dat$CAS)
    
    sample_dat <- filter(chemicalSummary, site %in% temp_site) %>%
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
           med_ndtect_new_rank = rank(-med_ndetect_new),
           max_ndetect_new_rank = rank(-max_ndetect_new)) %>%
    select(site, contains('rank'))
  
  nchems <- left_join(nchems, nchem_ranks)
  
  return(nchems)
}

get_conc_sites <- function(chemicalSummary_conc, all_samples) {
  site_conc <- chemicalSummary_conc %>%
    group_by(site, date) %>%
    summarize(sum_conc = sum(EAR), 
              mean_conc = mean(EAR))
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    filter()
  
  site_conc <- full_join(site_conc, samples) %>%
    mutate(sum_conc = ifelse(is.na(sum_conc), 0, sum_conc),
           mean_conc = ifelse(is.na(mean_conc), 0, mean_conc))
  
  site_conc_rank <- mutate(site_conc, 
                           med_sum_conc = median(sum_conc),
                           max_sum_conc = max(sum_conc),
                           med_mean_conc = median(mean_conc),
                           max_mean_conc = max(mean_conc)) %>%
    mutate(med_sum_conc_rank = rank(-med_sum_conc),
           max_sum_conc_rank = rank(-max_sum_conc),
           med_mean_conc_rank = rank(-med_mean_conc),
           max_mean_conc_rank = rank(-max_mean_conc)) %>%
    select(site, contains('rank'))
  
  site_conc <- left_join(site_conc, site_conc_rank)
                           
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
  
  samples <- all_samples %>%
    mutate(date = as.POSIXct(date)) %>%
    filter()
  
  site_ear <- full_join(site_ear, samples) %>%
    mutate(mean_sumEAR = ifelse(is.na(mean_sumEAR), 0, mean_sumEAR),
           max_sumEAR = ifelse(is.na(max_sumEAR), 0, max_sumEAR))
  
  site_ear_rank <- mutate(site_ear, 
                           med_mean_sumEAR = median(mean_sumEAR),
                           max_mean_sumEAR = max(mean_sumEAR),
                           med_max_sumEAR = median(max_sumEAR),
                           max_max_sumEAR = max(max_sumEAR)) %>%
    mutate(med_mean_sumEAR_rank = rank(-med_mean_sumEAR),
           max_mean_sumEAR_rank = rank(-max_mean_sumEAR),
           med_max_sumEAR_rank = rank(-med_max_sumEAR),
           max_max_sumEAR_rank = rank(-max_max_sumEAR)) %>%
    select(site, contains('rank'))
  
  site_ear <- left_join(site_ear, site_ear_rank)
  
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
           max_bench = ifelse(is.na(max_bench), 0, max_bench))
  
  site_bench_rank <- mutate(site_bench, 
                           med_mean_bench = median(mean_bench),
                           max_mean_bench = max(mean_bench),
                           med_max_bench= median(max_bench),
                           max_max_bench = max(max_bench)) %>%
    mutate(med_mean_bench_rank = rank(-med_mean_bench),
           max_mean_bench_rank = rank(-max_mean_bench),
           med_max_bench_rank = rank(-med_max_bench),
           max_max_bench_rank = rank(-max_max_bench)) %>%
    select(site, contains('rank'))
  
  site_bench <- left_join(site_bench, site_bench_rank)
  
  return(site_bench)
  
}
