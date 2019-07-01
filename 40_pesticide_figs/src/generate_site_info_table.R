# script to generate "sites" table for paper

gather_site_info <- function(target_name, sites, sum_conc, unique_chems) {
  site_conc_summary <- sum_conc %>%
    group_by(SiteID) %>%
    summarize(sample_count = n(),
              mean_dconc = round(mean(sum_conc_detect), 2),
              sd_dconc = round(sd(sum_conc_detect), 2),
              min_dconc = round(min(sum_conc_detect), 2),
              max_dconc = round(max(sum_conc_detect), 2), 
              mean_samp_nchem = round(mean(n_detected), 0), 
              min_samp_nchem = min(n_detected), 
              max_samp_nchem = max(n_detected),
              n_with_neonics = sum(neonic_meas))
  
  n_unique_chems <- unique_chems %>%
    group_by(SiteID) %>%
    summarize(n_unique_chems = sum(n_added_on_date))
  
  site_dat <- sites %>%
    select(SiteID = site_no, Site.name, shortName, lat = dec_lat_va, lon = dec_long_va, 
           perc_ag = `Ag..total`, perc_urban = Urban, dominant_lulc = `Dominant.land.use.`)
  
  site_table <- left_join(site_conc_summary, n_unique_chems)
  
  site_table <- left_join(site_dat, site_table) %>%
    select(`USGS Site Number` = SiteID, `Site Name` = Site.name, `Site Short Name` = shortName, Latitude = lat, Longitude = lon, 
           `% Agriculture` = perc_ag, `% Urban` = perc_urban, `Dominant Land Use/Land Cover` = dominant_lulc, `Sample Count` = sample_count,
           `Neonic Sample Count` = n_with_neonics)
  
  write.csv(site_table, target_name, row.names = F)
}

gather_chem_info <- function(out_file, chem_dls, chemicalSummary_conc, chem_info, parents) {
  
  # chemical table for supplement
  # includes chem name + CAS, n_detect_frequency (n_detect), median DL, median (min, max) concentration, 
  
  detect <- chem_dls %>%
    mutate(`Detection frequency, %` = round(100*(n_detect/(n_detect+n_bdl)), 0),
           `Median detection limit` = value) %>%
    select(pCode, `Detection frequency, %`, `Median detection limit`)
  
  obs <- chemicalSummary_conc %>%
    group_by(CAS) %>%
    summarise(median = round(median(EAR), 3),
              min = round(min(EAR), 3), 
              max = round(max(EAR), 3),
              `N sites` = length(unique(site)),
              `N months` = length(unique(lubridate::month(date)))) %>%
    mutate('Observed concentration, median (min, max)' = paste0(median, ' (', min, ', ', max, ')'))
  
  info <- left_join(detect, select(dataRetrieval::parameterCdFile, pCode = parameter_cd, CAS = casrn)) %>%
    left_join(chem_info) %>%
    left_join(select(parents, pCode, MlWt, `Parent compound` = parent_pesticide))
  
  missing_chems <- parameterCdFile[which(parameterCdFile$parameter_cd %in% info$pCode[is.na(info$`Chemical Name`)]),]
  missing_chem_names <- gsub(', water.*', '', missing_chems$parameter_nm) 
  
  info$`Chemical Name`[is.na(info$`Chemical Name`)] <- missing_chem_names
  
  info <- info %>%
    select(`Chemical Name`, CAS, `Molecular weight` = MlWt, Class, `Parent compound`, pCode) %>%
    mutate(`Parent compound` = ifelse(is.na(`Parent compound`), `Chemical Name`, `Parent compound`))
  
  
  out <- info %>%
    left_join(detect) %>%
    left_join(obs) %>%
    select(-pCode, -median, -min, -max) %>%
    arrange(-`Detection frequency, %`)
  
  write.csv(out, out_file, row.names = FALSE)
  
  
}
