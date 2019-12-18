clean_pesticides <- function(dat) {

  dat$source <- 'pesticides_s2437'
  
  # sum cis and trans permethrin, and update metadata
  no_permethrin <- filter(dat, !pCode %in% c('68708', '68769'))
  permethrin <- filter(dat, pCode %in% c('68708', '68769')) %>%
    group_by(pdate, SiteID, samp_type_cd, medium_cd, sample_dt, sample_tm) %>%
    summarize(value = sum(value),
              remark_cd = ifelse(all(is.na(remark_cd)), NA, '<')) %>%
    mutate(pCode = 'fake_1')
  
  dat <- bind_rows(no_permethrin, permethrin)
  
  return(dat)
}

clean_neonics <- function(dat) {
  
  dat$source <- 'neonic'
  return(dat)
}

clean_glyphosate <- function(dat) {
  dat$source <- 'glyphosate'
  return(dat)
}