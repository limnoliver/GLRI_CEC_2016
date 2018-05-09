library(dataRetrieval)
library(tidyr)
library(dplyr)

# This should include both schedule 4433 (wastewater indicators)
# and 2437 pesticides

get_NWIS <- function(tracking, schedule_pCodes, pCodesExclude){

  ## define sites for data retrieval ##
  sites <- zeroPad(unique(tracking$SiteID),8) #Use all sites in list for now

  df <- readNWISqw(siteNumbers = c("04157000",sites), parameterCd = "All",
                    startDate = '2015-10-01',
                    endDate = '2016-09-30')
  
  df_sub <- select(df,
                   pdate = startDateTime,
                   SiteID = site_no,
                   pCode = parm_cd,
                   value = result_va,
                   remark_cd, 
                   samp_type_cd, 
                   sample_dt,
                   sample_tm)
  #Saginaw!!!
  df_sub$SiteID[df_sub$SiteID == "04157005"] <- "04157000"
  
  schedule_pCodes <- filter(schedule_pCodes, 
                            !grepl("surr",schedule_pCodes$`Parameter Name`),
                            !(`Parameter Code` %in% pCodesExclude))
  
  pCodesToUse <- c(schedule_pCodes$`Parameter Code`,"99960")#, "62722", "62649")
  
  nwis_filtered <- filter(df_sub, pCode %in% pCodesToUse) %>%
    filter(samp_type_cd == "9") %>%
    select(-samp_type_cd)
  
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
  
  return(pCodeInfo)
}

