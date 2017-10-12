library(dplyr)
library(tidyr)
library(readr)
library(dataRetrieval)
library(openxlsx)
library(toxEval)

merged_NWIS <- function(tracking, NWIS, neonic, pCodeInfo, schedule_pCodes){

  tracking <- filter(tracking, SampleTypeCode == 9)
  
  just_neonic_data <- neonic %>%
    select( site = USGS.Site.ID, pdate, NWISRecordNumber,
           Acetamiprid, 
           Clothianidin,
           Dinotefuran,
           Imidacloprid,
           Thiacloprid,
           Thiamethoxam) %>%
    gather(chemical, value, -site, -pdate, -NWISRecordNumber) %>%
    mutate(site = zeroPad(site,8))
  
  just_neonic_remarks <- neonic %>%
    select(site = USGS.Site.ID, pdate, NWISRecordNumber,
           R_Acetamiprid, 
           R_Clothianidin,
           R_Dinotefuran,
           R_Imidacloprid,
           R_Thiacloprid,
           R_Thiamethoxam
    ) %>%
    gather(chemical_rk, remark_cd, -site, -pdate, -NWISRecordNumber) %>%
    mutate(chemical = gsub("R_","",chemical_rk)) %>%
    select(-chemical_rk)%>%
    mutate(site = zeroPad(site,8))

  just_neonic <- left_join(just_neonic_data, 
                           just_neonic_remarks, by=c("site","pdate","chemical","NWISRecordNumber"))

  just_NWIS <- select(NWIS, site=SiteID, NWISRecordNumber, pdate, pCode, value) %>%
    left_join(select(pCodeInfo, pCode=parameter_cd, chemical=casrn), by="pCode") %>%
    filter(pCode %in% schedule_pCodes$`Parameter Code`) %>%
    filter(NWISRecordNumber %in% tracking$NWISRecordNumber)
  
    # right_join(select(tracking, SiteID, pdate), by=c("site"="SiteID","pdate"))
  just_NWIS$site[just_NWIS$site == "04157005"] <- "04157000"  
  
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
                                  "Thiamethoxam") ,
                          Class = "Neonic",stringsAsFactors = FALSE)
  return(special_cas)
}

create_chemData <- function(neonic_NWIS, special_cas, pCodeInfo){

  chem_data <- neonic_NWIS %>%
    select(SiteID = site,
           `Sample Date` = pdate,
           Value = value,remark_cd,
           CAS = chemical,
           pCode) %>%
    filter(!is.na(Value),
           !is.na(CAS)) %>%
    left_join(special_cas, by="CAS") %>%
    left_join(select(pCodeInfo, pCode=parameter_cd, units=parameter_units), by="pCode") %>%
    mutate(Value = Value/1000)
    
  # chem_data$Value[is.na(chem_data$units)] <- chem_data$Value[is.na(chem_data$units)]/1000  # Neonics
  # chem_data$Value[chem_data$units == "ng/l"] <- chem_data$Value[chem_data$units == "ng/l"]/1000

  chem_data$CAS[!is.na(chem_data$casrn)] <- chem_data$casrn[!is.na(chem_data$casrn)]
  
  chem_data <- select(chem_data, -casrn, -Class) %>%
    distinct()

  return(chem_data)
  
  
}


create_tox_chemInfo <- function(neonic_NWIS, special_cas, pCodeInfo, classes){

  chem_data <- neonic_NWIS %>%
    select(SiteID = site,
           `Sample Date` = pdate,
           Value = value,remark_cd,
           CAS = chemical,
           pCode) %>%
    filter(!is.na(Value),
           !is.na(CAS)) %>%
    left_join(special_cas, by="CAS") 
  
  chem_data$CAS[!is.na(chem_data$casrn)] <- chem_data$casrn[!is.na(chem_data$casrn)]
  
  chem_info <- select(chem_data, CAS, pCode) %>%
    distinct() %>%
    left_join(select(pCodeInfo, CAS = casrn,
                     `Chemical Name` = parameter_nm,
                     pCode = parameter_cd), by=c("CAS","pCode")) %>%
    left_join(special_cas, by=c("CAS"="casrn"))
  
  chem_info$`Chemical Name`[!is.na(chem_info$CAS.y)] <- chem_info$CAS.y[!is.na(chem_info$CAS.y)]
  
  chem_info <- select(chem_info, -CAS.y) %>%
    left_join(select(classes, CAS, class1 = Class), by="CAS")

  chem_info$Class[is.na(chem_info$Class)] <- chem_info$class1[is.na(chem_info$Class)]
  
  chem_info <- select(chem_info, -class1, -pCode)
  
  return(chem_info)
}

create_tox_siteInfo <- function(sites){

  siteInfo <- sites %>%
    select(SiteID = site_no,
           `Short Name` = shortName,
           dec_lat = dec_lat_va,
           dec_lon = dec_long_va,
           map_nm,
           station_nm)%>%
    mutate(site_grouping = "All")
  
  siteInfo$`Short Name`[is.na(siteInfo$`Short Name`)] <- siteInfo$map_nm[is.na(siteInfo$`Short Name`)]
  siteInfo <- select(siteInfo, -map_nm)
  
  siteInfo$`Short Name`[is.na(siteInfo$`Short Name`)] <- siteInfo$station_nm[is.na(siteInfo$`Short Name`)]
  siteInfo <- select(siteInfo, -station_nm)
  
  return(siteInfo)
  
}
  
create_toxExcel <- function(chem_data, chem_info, site_info, file_out){
  

  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)

}

get_chem_sum <- function(chem_info, chem_data, site_info){
  ACClong <- get_ACC(chem_info$CAS)
  ACClong <- remove_flags(ACClong)
  
  cleaned_ep <- clean_endPoint_info(endPointInfo)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(ACClong,
                                          filtered_ep,
                                          chem_data, 
                                          site_info, 
                                          chem_info)
  return(chemicalSummary)
}
