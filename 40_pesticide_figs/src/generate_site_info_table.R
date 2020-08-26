# script to generate "sites" table for paper

gather_site_info <- function(target_name, sites, neonics) {
  # site_conc_summary <- sum_conc %>%
  #   group_by(SiteID) %>%
  #   summarize(sample_count = n(),
  #             mean_dconc = round(mean(sum_conc_detect), 2),
  #             sd_dconc = round(sd(sum_conc_detect), 2),
  #             min_dconc = round(min(sum_conc_detect), 2),
  #             max_dconc = round(max(sum_conc_detect), 2), 
  #             mean_samp_nchem = round(mean(n_detected), 0), 
  #             min_samp_nchem = min(n_detected), 
  #             max_samp_nchem = max(n_detected),
  #             n_with_neonics = sum(neonic_meas))
  # 
  # n_unique_chems <- unique_chems %>%
  #   group_by(SiteID) %>%
  #   summarize(n_unique_chems = sum(n_added_on_date))

  neonic_count <- neonics %>%
    filter(!is.na(value)) %>%
    group_by(SiteID) %>%
    summarize(`Neonic sample count` = length(unique(sample_dt)))
  

  site_dat <- sites %>%
    select(SiteID = site_no, Site.name, shortName, lat = dec_lat_va, lon = dec_long_va, 
           perc_crop = `Ag..crops`, perc_pasture = `Ag..pasture..hay.`, perc_urban = Urban, dominant_lulc = `Dominant.land.use.`)
  
  #site_table <- left_join(site_conc_summary, n_unique_chems)
  site_table <- sites %>%
    left_join(neonic_count, by = c('site_no' = 'SiteID')) %>% 
    select(`USGS Site Number` = site_no, `Site Name` = Site.name, `Site Short Name` = shortName, 
           Latitude = dec_lat_va, Longitude = dec_long_va, 
           `% Crops` = `Ag..crops`, `% Pasture` = `Ag..pasture..hay.`, `% Urban` = Urban, 
           `% Forest` = Forest, `% Wetland` = `Water..wetland`,
           `Dominant Land Use/Land Cover` = `Dominant.land.use.`,
           `Neonic sample count`) %>%
    mutate(`Neonic sample count` = ifelse(is.na(`Neonic sample count`),0, `Neonic sample count`)) %>%
    mutate(`Sample count` = ifelse(`Site Short Name` %in% 'Maumee', 18, 12))
  
  write.csv(site_table, target_name, row.names = F)
}

gather_chem_info <- function(out_file, chem_dls, reduced_dat, chem_crosswalk, chem_info_all) {
  
  # chemical table for supplement
  # includes chem name + CAS, n_detect_frequency (n_detect), median DL, median (min, max) concentration, 
  
  detect <- chem_dls %>%
    mutate(`Detection frequency` = paste0(round(100*(n_detect/(n_detect+n_bdl)), 0), ' (', n_detect, '/', n_all, ')'),
           `Detection limit` = paste0(value, ' (', min, ', ', max, ')'),
           ) %>%
    select(pCode, `Detection frequency`, `Detection limit`)
  
  obs <- reduced_dat %>%
    group_by(chnm, CAS, Class, parent_pesticide) %>%
    summarise(Concentration = if_else(length(which(EAR>0)) == 0, 'NA', paste0(round(median(EAR[EAR>0]), 3), ' (', round(max(EAR), 3), ')')),
              Sites = paste0(length(unique(site[EAR>0])), '/' , length(unique(site))),
              Months = paste0(length(unique(lubridate::month(date[EAR>0]))), '/', length(unique(lubridate::month(date))))) %>%
    left_join(select(chem_crosswalk, CAS, pCode))
  
  out <- left_join(detect, chem_crosswalk) %>%
    left_join(obs) %>%
    left_join(select(chem_info_all, CAS, pCode))
  
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
