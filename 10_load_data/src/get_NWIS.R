library(dataRetrieval)
library(tidyr)
library(dplyr)

# This should include both schedule 4433 (wastewater indicators)
# and 2437 pesticides

get_NWIS <- function(tracking){

  ## define sites for data retrieval ##
  sites <- zeroPad(unique(tracking$SiteID),8) #Use all sites in list for now

  df <- readNWISqw(siteNumbers = c("04157005",sites), parameterCd = "All",
                    startDate = '2015-10-01',
                    endDate = '2016-09-30')
  
  df_sub <- select(df,
                   pdate = startDateTime,
                   SiteID = site_no,
                   pCode = parm_cd,
                   value = result_va,
                   remark_cd, 
                   samp_type_cd,
                   medium_cd,
                   sample_dt,
                   sample_tm)
  #Saginaw!!!
  df_sub$SiteID[df_sub$SiteID == "04157005"] <- "04157000"
  
  #df_sub <- filter(df_sub, samp_type_cd == '9') %>%
  #  select(-samp_type_cd)
  
  return(df_sub)
  
}

filter_pesticides <- function(NWIS, schedule_pCodes, pCodesExclude) {

  schedule_pCodes <- filter(schedule_pCodes, 
                            schedule %in% 2437, 
                            !grepl("surr",schedule_pCodes$`Parameter Name`),
                            !(`Parameter Code` %in% pCodesExclude), 
                            !(`Parameter Code` %in% '68426')) # exclude imidacloprid for now, add back in once have unique date/times
  
  pCodesToUse <- c(schedule_pCodes$`Parameter Code`)

  nwis_filtered <- filter(NWIS, pCode %in% pCodesToUse)
  
  pest_sum <- select(nwis_filtered, SiteID, pdate) %>%
    distinct()
  
  imidacloprid <- filter(NWIS, pCode %in% '68426')
  imidacloprid <- left_join(pest_sum, imidacloprid, by = c('SiteID', 'pdate'))
  
  all_dat <- bind_rows(nwis_filtered, imidacloprid)
  
  return(all_dat)
  
}


filter_neonics <- function(NWIS) {
  
  pCodes_neonics <- c('68221', '68245', '68302', '68379', '68485') # adds the neonics except Imidacloprid which is also in schedule
  # reduce to neonics
  
  neonics <- filter(NWIS, pCode %in% pCodes_neonics)
  
  # find unique site/date combinations to pull Imidocloprid from neonic measurements (not schedule)
  neonic_sum <- select(neonics, SiteID, pdate) %>%
    distinct()
  
  imidacloprid <- filter(NWIS, pCode %in% '68426') 
  imidacloprid <- left_join(neonic_sum, imidacloprid, by = c('SiteID', 'pdate'))
    
  
  neonics_all <- bind_rows(neonics, imidacloprid)
  
  return(neonics_all)
  site_sum <- group_by(neonic_sum, SiteID) %>%
    summarize(n = n())
  
  test_missing <- filter(neonics, SiteID == "04208000") %>%
    group_by(pdate) %>%
    summarize(n = n())
  
  test_site <- filter(NWIS, SiteID == "04208000") %>%
    group_by(pdate) %>%
    summarize(n = n())
  
  test_datetime <- filter(NWIS, pdate == as.POSIXct('2015-10-01 15:15:00'))

}

filter_owc <- function(NWIS, schedule_pCodes) {
  
  schedule_pCodes <- filter(schedule_pCodes, 
                            schedule %in% 4433, 
                            !grepl("surr",schedule_pCodes$`Parameter Name`))
  
  pCodesToUse <- c(schedule_pCodes$`Parameter Code`)
  
  nwis_filtered <- filter(NWIS, pCode %in% pCodesToUse)
  
  return(nwis_filtered)
}

filter_glyphosate <- function(NWIS) {
  pCodes_glyphosate <- c("62722", "62649","99960") #codes for glyphosate and degradate 

  # reduce to glyphosate
  glyphosate <- filter(NWIS, pCode %in% pCodes_glyphosate)
  
  return(glyphosate)
  
}

merge_tracking <- function() {
  nwis_tracking <- left_join(nwis_filtered, 
                             tracking, 
                             by=c("sample_dt"="Date",
                                  "sample_tm"="Time",
                                  "SiteID")) %>%
    rename(pdate = pdate.x)
  
  
  # x <- tracking[which(!(tracking$NWISRecordNumber %in% unique(nwis_tracking$NWISRecordNumber))),]
  # y <- filter(nwis_tracking, is.na(NWISRecordNumber))
  ##########################
  
  return(nwis_tracking)
  
}


get_pCode_exclude <- function(path_to_exclude){

  pCodes <- read.csv(path_to_exclude, header = FALSE,
                     stringsAsFactors = FALSE)
  return(pCodes$V1)
}

get_pCode <- function(NWIS){

  pCodes <- unique(NWIS$pCode)
  pCodeInfo <- readNWISpCode(pCodes) 
  pCodeInfo$casrn[pCodeInfo$parameter_cd %in% '68574'] <- '56611-55-3'
  
  
  
  return(pCodeInfo)
}

