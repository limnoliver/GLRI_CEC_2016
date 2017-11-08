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
  
  # Making things right:
  tracking$SiteID[tracking$SiteID == "04157005"] <- "04157000"
  tracking <- filter(tracking, SampleTypeCode == "9")
  tracking$Date[tracking$NWISRecordNumber == "01601472"] <- as.Date("2016-04-07")
  tracking$Date[tracking$NWISRecordNumber == "01600121"] <- as.Date("2016-02-02")
  tracking$SiteID[tracking$SiteID == "4249000"] <- "04249000"
  tracking$Time[tracking$NWISRecordNumber == "01600388"] <- "10:20"
  tracking$Time[tracking$NWISRecordNumber == "01600462"] <- "09:10"
  tracking$Date[tracking$NWISRecordNumber == "01600384"] <- as.Date("2015-12-09")
  tracking$Date[tracking$NWISRecordNumber == "01600381"] <- as.Date("2015-12-08")
  tracking$Time[tracking$NWISRecordNumber == "01601976"] <-  "08:00"
  tracking$Date[tracking$NWISRecordNumber == "01600102"] <- as.Date("2015-11-12")
  tracking$Time[tracking$NWISRecordNumber == "01600461"] <- "12:10"
  tracking$Time[tracking$NWISRecordNumber == "016002001"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01600207"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01600236"] <- "12:00"
  tracking$Time[tracking$NWISRecordNumber == "01600484"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01600973"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01601028"] <- "12:00"
  tracking$Date[tracking$NWISRecordNumber == "01601028"] <- as.Date("2016-05-10")
  tracking$Time[tracking$NWISRecordNumber == "01601556"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01601641"] <- "10:00"
  tracking$Time[tracking$NWISRecordNumber == "01600618"] <- "13:00"
  tracking$Time[tracking$NWISRecordNumber == "01600700"] <- "10:00"
  
  schedule_pCodes <- filter(schedule_pCodes, 
                            !grepl("surr",schedule_pCodes$`Parameter Name`),
                            !(`Parameter Code` %in% pCodesExclude))
  
  nwis_filtered <- filter(df_sub, pCode %in% schedule_pCodes$`Parameter Code`) %>%
    filter(samp_type_cd == "9") %>%
    select(-samp_type_cd)
  
  nwis_tracking <- left_join(nwis_filtered, 
                             tracking, 
                             by=c("sample_dt"="Date",
                                  "sample_tm"="Time",
                                  "SiteID")) %>%
    rename(pdate = pdate.x)
  
  x <- tracking[which(!(tracking$NWISRecordNumber %in% unique(nwis_tracking$NWISRecordNumber))),]
  y <- filter(nwis_tracking, is.na(NWISRecordNumber))
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

get_AOPs <- function(){
  path_to_tox <-  system.file("extdata", package="toxEval")
  file_name <- "AOP_crosswalk.csv"
  full_path <- file.path(path_to_tox, file_name)
  
  AOP_crosswalk <- read.csv(full_path, stringsAsFactors = FALSE)
  
  AOP_crosswalk <- select(AOP_crosswalk, 
                          gene_symbol=Target.Gene.Symbol, 
                          AOP=AOP.name)
  
  AOP <- data.frame(gene_symbol = unique(AOP_crosswalk$gene_symbol),
                    AOP = "",
                    stringsAsFactors = FALSE)
  for(gene in AOP$gene_symbol){
    AOP$AOP[AOP$gene_symbol %in% gene] <- paste(AOP_crosswalk$AOP[AOP_crosswalk$gene_symbol %in% gene],collapse = ", ")
  }
  
  return(AOP)
}