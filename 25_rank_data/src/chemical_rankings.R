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