library(dplyr)
library(USGSHydroTools)
library(tidyr)
library(readr)
library(dataRetrieval)
library(openxlsx)

merged_neonic_flow <- function(neonic, flow, sites){

  flow$Q <- ifelse(!is.na(flow$X_00060_00000),flow$X_00060_00000,flow$X_00060_00003)
  
  listQ <- list()
  
  for(site in sites$USGS.station.number){
    
    if(site == "04157005") site <- "04157000"
    
    sub_data <- filter(neonic, USGS.station.number == site)
    dfQ2 <- filter(flow, site_no == site)
    
    if(nrow(dfQ2) > 0){
      subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = sub_data, starttime = "pdate",
                             times=6,units="hrs",stats.return = c("mean","max","min","difference"),
                             subdatesvar = "",out.varname = "Q")
      subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                             times=0,units="hrs",stats.return = c("nearprev"),
                             subdatesvar = "",out.varname = "Q")
      subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                             times=1,units="hrs",stats.return = c("mean"),
                             subdatesvar = "",out.varname = "Q")
      listQ[[site]] <- subdfNeonic
    }
    
  }

  dfNeonic2 <- listQ[[1]]

  for(i in 2:length(listQ)){
    dfNeonic2 <- bind_rows(dfNeonic2,listQ[[i]])
  }
  
  return(dfNeonic2)
}

merged_NWIS <- function(tracking, NWIS, neonic_flow, pCodeInfo){
  
  just_neonic_data <- neonic_flow %>%
    select(Sample, site = USGS.station.number, pdate, 
           Acetamiprid, 
           Clothianidin,
           Dinotefuran,
           Imidacloprid,
           Thiacloprid,
           Thiamethoxam) %>%
    gather(chemical, value, -Sample, -site, -pdate)
  
  just_neonic_remarks <- neonic_flow %>%
    select(Sample, site = USGS.station.number, pdate, 
           R_Acetamiprid, 
           R_Clothianidin,
           R_Dinotefuran,
           R_Imidacloprid,
           R_Thiacloprid,
           R_Thiamethoxam
    ) %>%
    gather(chemical_rk, remark_cd, -Sample, -site, -pdate) %>%
    mutate(chemical = gsub("R_","",chemical_rk)) %>%
    select(-chemical_rk)

  just_neonic <- left_join(just_neonic_data, 
                           just_neonic_remarks, by=c("Sample","site","pdate","chemical"))
  
  just_NWIS <- select(NWIS, site=SiteID, Sample = NWISRecordNumber, pdate, pCode, value) %>%
    left_join(select(pCodeInfo, pCode=parameter_cd, chemical=casrn), by="pCode") 
  
  nwis_neonic <- bind_rows(just_neonic, just_NWIS)
  
  return(nwis_neonic)
}

get_special_cas <- function(){
  special_cas <- data.frame(casrn = c("135410-20-7",
                                    "210880-92-5",
                                    "165252-70-0",
                                    "138261-41-3",
                                    "111988-49-9",
                                    "153719-23-4"),
                          CAS = c("Acetamiprid",
                                  "Clothianidin",
                                  "Dinotefuran",
                                  "Imidacloprid",
                                  "Thiacloprid",
                                  "Thiamethoxam") ,stringsAsFactors = FALSE)
  return(special_cas)
}

create_chemData <- function(neonic_NWIS, special_cas){

  chem_data <- neonic_NWIS %>%
    select(SiteID = site,
           `Sample Date` = pdate,
           Value = value,remark_cd,
           CAS = chemical) %>%
    filter(!is.na(Value),
           !is.na(CAS)) %>%
    left_join(special_cas, by="CAS") %>%
    mutate(Value = Value*1000) #ng -> ug
  
  chem_data$CAS[!is.na(chem_data$casrn)] <- chem_data$casrn[!is.na(chem_data$casrn)]
  
  chem_data <- select(chem_data, -casrn)

  return(chem_data)
  
  
}


create_tox_chemInfo <- function(chem_data, special_cas, pCodeInfo){

  chem_info <- select(chem_data, CAS) %>%
    distinct() %>%
    left_join(select(pCodeInfo, CAS = casrn,
                     `Chemical Name` = parameter_nm), by="CAS") %>%
    left_join(special_cas, by=c("CAS"="casrn"))
  
  chem_info$`Chemical Name`[!is.na(chem_info$CAS.y)] <- chem_info$CAS.y[!is.na(chem_info$CAS.y)]
  
  chem_info <- select(chem_info, -CAS.y) %>%
    mutate(Class = "Neonic")

  return(chem_info)
}

create_tox_siteInfo <- function(sites){

  siteInfo <- sites %>%
    select(SiteID = USGS.station.number,
           `Short Name` = shortName)
  
  USGS_site <- readNWISsite(siteInfo$SiteID)
  
  siteInfo <- siteInfo %>%
    left_join(select(USGS_site, SiteID = site_no, 
                     dec_lat = dec_lat_va,
                     dec_lon = dec_long_va), by="SiteID") %>%
    mutate(site_grouping = "All")
  
  return(siteInfo)
  
}
  
create_toxExcel <- function(chem_data, chem_info, site_info, file_out){
  

  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)

}
