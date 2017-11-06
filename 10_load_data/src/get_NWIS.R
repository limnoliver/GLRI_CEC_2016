library(dataRetrieval)
library(tidyr)
library(dplyr)

# This should include both schedule 4433 (wastewater indicators)
# and 2437 pesticides

get_NWIS <- function(tracking){

  ## define sites for data retrieval ##
  sites <- zeroPad(unique(tracking$SiteID),8) #Use all sites in list for now

  df <- readWQPdata(siteNumbers = c("USGS-04157000",paste0("USGS-",sites)), 
                      startDate = '2015-10-01',
                      endDate = '2016-09-30')

  df$remark_cd <- ""
  df$remark_cd[grep("^Not Detected$",df$ResultDetectionConditionText)] <- "<"
  df$remark_cd[grep("Detected Not Quantified",df$ResultDetectionConditionText)] <- "<"
  
  
  df$occur <-  df$remark_cd != "<"
  df$occur <- ifelse(is.na(df$occur),TRUE,df$occur)
  
  df_sub <- df %>%
    select(pdate = ActivityStartDateTime, FullsiteID = MonitoringLocationIdentifier, ActivityIdentifier,
           CharacteristicName,pCode=USGSPCode,value = ResultMeasureValue,remark_cd,units = ResultMeasure.MeasureUnitCode,
           ActivityTypeCode, HydrologicCondition, HydrologicEvent,ResultValueTypeName,
           lab.comments = ResultLaboratoryCommentText,
           detection.type = DetectionQuantitationLimitTypeName,
           detection.limit = DetectionQuantitationLimitMeasure.MeasureValue,
           detection.units = DetectionQuantitationLimitMeasure.MeasureUnitCode) %>%
    mutate(SiteID = gsub(pattern = "USGS-","",FullsiteID),
           NWISRecordNumber = sapply(strsplit(ActivityIdentifier,"\\."), function(x) x[3])) %>%
    select(-ActivityIdentifier, -FullsiteID)

  #Saginaw!!!
  df_sub$SiteID[df_sub$SiteID == "04157005"] <- "04157000"

  return(df_sub)
  
}

get_pCode_exclude <- function(path_to_exclude){

  pCodes <- read.csv(path_to_exclude, header = FALSE,
                     stringsAsFactors = FALSE)
  return(pCodes$V1)
}

get_pCode <- function(NWIS, pCodesExclude){

  pCodes <- unique(NWIS$pCode)
  pCodeInfo <- readNWISpCode(pCodes) %>%
    filter(!(parameter_cd %in% pCodesExclude))
  
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