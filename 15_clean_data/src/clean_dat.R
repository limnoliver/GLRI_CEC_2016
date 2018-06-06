clean_pesticides <- function(dat) {
  # Milwaukee River has tracking date of 2016-04-07 but sampling date
  # of 2016-04-05. Will change in data so as not to affect any other data that is merged
  # with tracking
  
  dat$sample_dt[dat$SiteID %in% '04087170' & dat$sample_dt %in% as.Date('2016-04-05')] <- as.Date('2016-04-07')
  dat$source <- 'pesticides_s2437'
  return(dat)
}

clean_neonics <- function(dat) {
  
  # Maumee river site has tracking date of 2016-05-10 but sampling date of 2016-05-03. 
  # Hladik has 2016-05-03 published as sample date, but this didn't seem to trip up pesticides or glyphosate data,
  # so am just changing neonic date. 
  dat$sample_dt[dat$SiteID %in% '04193500' & dat$sample_dt %in% as.Date('2016-05-03')] <- as.Date('2016-05-10')
  
  dat$source <- 'neonic'
  return(dat)
}

clean_glyphosate <- function(dat) {
  dat$sample_dt[dat$SiteID %in% '04087170' & dat$sample_dt %in% as.Date('2016-04-05')] <- as.Date('2016-04-07')
  dat$source <- 'glyphosate'
  return(dat)
}