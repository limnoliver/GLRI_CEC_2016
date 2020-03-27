clean_pesticides <- function(dat) {

  dat$source <- 'pesticides_s2437'
  
  # sum cis and trans permethrin, and update metadata
  no_permethrin <- filter(dat, !pCode %in% c('68708', '68769'))
  permethrin <- filter(dat, pCode %in% c('68708', '68769')) %>%
    group_by(pdate, SiteID, samp_type_cd, medium_cd, sample_dt, sample_tm) %>%
    summarize(value = sum(value),
              remark_cd = ifelse(all(is.na(remark_cd)), NA, '<')) %>%
    mutate(pCode = 'fake_1')
  
  dat <- bind_rows(ungroup(no_permethrin), ungroup(permethrin))
  
  return(dat)
}

clean_neonics <- function(dat) {
  
  dat$source <- 'neonic'
  return(dat)
}

clean_glyphosate <- function(dat) {
  dat$source <- 'glyphosate'
  
  # drop the full glyphosate analysis + degradate, use the immunoassay
  # eventually should estimate degradate based on glyphosate ~ degradate relationship - see Mahler's work on topic
  glyphosate_drop <- c("62722", "62649")
  
  # find ratio between degradate (62649) and parent
  immuno <- filter(dat, pCode %in% glyphosate_drop) %>%
    select(pdate, SiteID, pCode, value, remark_cd) %>%
    filter(!remark_cd %in% '<') %>%
    select(-remark_cd) %>%
    tidyr::pivot_wider(id_cols = c(pdate, SiteID, pCode), names_from = pCode, values_from = value) %>%
    filter(!is.na(`62649`)) %>% filter(!is.na(`62722`))
  
  # ratio of deg to parent for multiplier
  ratio <- median(immuno$`62649`/immuno$`62722`)
  
  fixed_dat <- filter(dat, !(pCode %in% glyphosate_drop))
  # keeping in "<" values so that the appropriate calculations for proportion detected
  # are kept. Do not adjust values, though. Remove this compound from DL analysis
  deg_fixed <- fixed_dat %>%
    mutate(pCode = '62649',
           value = ifelse(remark_cd %in% '<', value, value*ratio))
  
  fixed_dat <- bind_rows(fixed_dat, deg_fixed)
    
  return(fixed_dat)
}