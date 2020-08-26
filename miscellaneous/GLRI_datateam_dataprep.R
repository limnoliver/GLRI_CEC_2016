library(dplyr)

owc <- make('owc', remake_file = '10_load_data.yml') %>%
  select(sample_date = sample_dt, site_no = SiteID, pcode = pCode, value, remark_cd) %>%
  mutate(site_no = paste0('USGS-', site_no))

owc_pcodes <- filter(dataRetrieval::parameterCdFile, parameter_cd %in% unique(owc$pcode)) %>%
  rename(pcode = parameter_cd)

owc_sites <- make('sites', remake_file = '10_load_data.yml') %>%
  select(site_no, Site.name, shortName, dec_lat_va, dec_long_va)%>%
  mutate(site_no = paste0('USGS-', site_no)) %>%
  rename(site_name = Site.name, short_name = shortName, latitude = dec_lat_va, longitude = dec_long_va) %>%
  filter(site_no %in% unique(owc$site_no))
  

write.csv(owc, 'usgs_owc_2016.csv', row.names = FALSE)
write.csv(owc_pcodes, 'usgs_owc_pcodes_2016.csv', row.names = FALSE)
write.csv(owc_sites, 'usgs_owc_sites_2016.csv', row.names = FALSE)


pesticide <- make('pesticides', remake_file = '10_load_data.yml') %>%
  filter(!is.na(value)) %>%
  select(sample_date = sample_dt, site_no = SiteID, pcode = pCode, value, remark_cd) %>%
  mutate(site_no = paste0('USGS-', site_no))

pesticide_pcodes <- filter(dataRetrieval::parameterCdFile, parameter_cd %in% unique(pesticide$pcode)) %>%
  rename(pcode = parameter_cd) %>%
  mutate(pcode = paste0(' ', pcode))
 
pesticide_sites <- make('sites', remake_file = '10_load_data.yml') %>%
  select(site_no, Site.name, shortName, dec_lat_va, dec_long_va)%>%
  mutate(site_no = paste0('USGS-', site_no)) %>%
  rename(site_name = Site.name, short_name = shortName, 
         latitude = dec_lat_va, longitude = dec_long_va) %>%
  filter(site_no %in% unique(pesticide$site_no))

write.csv(pesticide, 'usgs_pesticide_2016.csv', row.names = FALSE)
write.csv(pesticide_pcodes, 'usgs_pesticide_pcodes_2016.csv', row.names = FALSE)
write.csv(pesticide_sites, 'usgs_pesticide_sites_2016.csv', row.names = FALSE)


