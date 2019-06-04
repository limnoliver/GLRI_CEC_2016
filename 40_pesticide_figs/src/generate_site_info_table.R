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