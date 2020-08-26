# find unique dates by site
# should be 12 per site

dat <- make('merged_dat')

head(dat)

dat.ordered <- arrange(dat, site, date, chnm)

dat.date.miss <- filter(dat.ordered, is.na(date))

# these appear to be neonics
# find all dates for bad river site for neonic measurements
bad.neonic.dates <- filter(dat, site == '04027000' & chnm == 'Thiacloprid')
bad.neonic.dates <- filter(dat, site == '04027000' & chnm == 'Acetamiprid')
bad.neonic.dates <- filter(dat, site == '04027000' & chnm == 'Acetamiprid')

# group by site and date and sum 

sum_dat <- dat %>%
  group_by(site, date) %>%
  summarize


head(dat.ordered)

# neonic vs tracking # 
# test new neonic pull from NWIS against tracking
tracking <- make('tracking')
tracking_sub <- filter(tracking, MediumCode == 'WS')
neonic <- make('neonic')
neonic_summary <- group_by(neonic, SiteID, sample_dt) %>%
  summarize(n = n())

neonic_j <- select(tracking_sub, Site, SiteID, Date, Neonics) %>%
  left_join(neonic_summary, by = c('SiteID', 'Date' = 'sample_dt')) %>%
  mutate(neonic_sample = ifelse(nchar(Neonics) %in% 6, Neonics, NA))

neonic_rj <- left_join(neonic_summary, select(tracking_sub, Site, SiteID, Date, Neonics), by = c('SiteID', 'sample_dt' = 'Date'))

neonic_track_no_nwis <- filter(neonic_j, !is.na(neonic_sample) & is.na(n))
neonic_nwis_no_track <- filter(neonic_rj, !is.na(n) & is.na(Neonics))  
head(neonic_j)

# find missing Indiana harbor canal neonic data site 04092750
ihc <- filter(NWIS, SiteID %in% '04092750')
ihc <- left_join(ihc, select(parameterCdFile, parameter_cd, parameter_nm), by = c('pCode' = 'parameter_cd'))

ihc[grep('thiacloprid', ihc$parameter_nm, ignore.case = T), ]

# pesticides vs tracking # 
pesticides <- make('pesticides')

pesticide_summary <- group_by(pesticides, SiteID, sample_dt) %>%
  summarize(n = n())

pesticide_j <- select(tracking_sub, Site, SiteID, Date, Comments) %>%
  left_join(pesticide_summary, by = c('SiteID', 'Date' = 'sample_dt'))

pesticide_rj <- left_join(pesticide_summary, select(tracking_sub, Site, SiteID, Date, NWISRecordNumber), by = c('SiteID', 'sample_dt' = 'Date'))

# look at pesticide data in NWIS that is not in tracking
pesticide_nwis_no_track <- filter(pesticide_rj, is.na(NWISRecordNumber))
pesticide_track_no_nwis <- filter(pesticide_j, is.na(n))
## glyphosates vs tracking ##

# find all glyphosate pcodes
NWIS <- make('NWIS')
all_pcodes <- parameterCdFile[parameterCdFile$parameter_cd %in% unique(NWIS$pCode), ]
glyphosate_pcodes <- all_pcodes[grep('glyphosate', all_pcodes$parameter_nm, ignore.case = T), ]

glyphosate <- make('glyphosate')

glyphosate_summary <- group_by(glyphosate, SiteID, sample_dt) %>%
  summarize(n = n())

glyphosate_j <- select(tracking_sub, Site, SiteID, Date, Glyphosate, GlyphosateComments) %>%
  left_join(glyphosate_summary, by = c('SiteID', 'Date' = 'sample_dt'))

glyphosate_rj <- left_join(glyphosate_summary, select(tracking_sub, Site, SiteID, Date, Glyphosate), by = c('SiteID', 'sample_dt' = 'Date'))

glyphosate_tracking_no_nwis <- filter(glyphosate_j, !is.na(Glyphosate) & is.na(n))
glyphosate_nwis_no_tracking <- filter(glyphosate_rj, is.na(Glyphosate))
