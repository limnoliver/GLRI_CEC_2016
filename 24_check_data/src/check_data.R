get_missing_chems <- function(conc_dat, special_cas) {
  
  # exclude chems without CAS numbers, also exclude neonics which aren't
  # supposed to be in every sample
  conc_dat <- conc_dat %>%
    filter(!(CAS %in% special_cas$casrn))
    
  
  n_by_sample <- group_by(conc_dat, site, date) %>%
    summarize(count = n())
  
  # max n is 215, 209 without the neonics
  n_max <- max(n_by_sample$count)
  
  # which samples don't have all 209?
  n_unique <- unique(conc_dat$CAS)
  
  missing <- group_by(conc_dat, SiteID, sample_dt) %>%
    summarize(names_chems = paste(CAS, collapse = ', '),
              n_missing_chems = length(which(!(n_unique %in% CAS))),
              cas_missing_chems = paste(n_unique[which(!(n_unique %in% CAS))], collapse = ', '))
  
  missing_unique_strings <- paste(unique(missing$cas_missing_chems), collapse = ", ")
  missing_unique_strings <- unlist(strsplit(missing_unique_strings, split=", "))
  missing_unique_strings <- unique(missing_unique_strings)
  missing_unique_strings <- missing_unique_strings[-which(missing_unique_strings %in% "")]
  
  n_per_chem <- conc_dat %>%
    group_by(CAS) %>%
    summarize(count = n(),
              count_det = length(which(!(remark_cd %in% "<"))), 
              det_rate = count_det/count)
  
  return(n_per_chem)
}

assess_neonics <- function(reduced_dat, special_cas) {
  neonics <- left_join(reduced_dat, select(parameterCdFile,pCode = parameter_cd, CAS = casrn)) %>%
    filter(CAS %in% special_cas$casrn) %>%
    group_by(CAS) %>%
    summarize(count = n(),
              count_det = length(which(!(remark_cd %in% "<"))), 
              det_rate = count_det/count) %>%
    left_join(select(special_cas, CAS = casrn, chnm = CAS))
  
  return(neonics)
  
  
}

assess_missing_chems <- function(missing_chems, EAR_dat, chem_data) {
  missing_chems <- filter(missing_chems, count < 198)
  
  missing_ears <- filter(EAR_dat, CAS %in% missing_chems$CAS) %>%
    group_by(CAS) %>%
    summarize(n_hits = sum(EAR>0.001),
              n_hit_rate = sum(EAR>0.001)/n())
  
  missing <- left_join(missing_chems, missing_ears, by = 'CAS')
    
    
}