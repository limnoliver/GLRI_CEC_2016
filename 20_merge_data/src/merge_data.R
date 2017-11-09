library(dplyr)
library(tidyr)
library(readr)
library(dataRetrieval)
library(openxlsx)
library(toxEval)

merged_NWIS <- function(tracking, NWIS, neonic, pCodeInfo){

  tracking <- filter(tracking, SampleTypeCode == 9)

  just_neonic_data <- neonic %>%
    select(site = USGS.Site.ID, NWISRecordNumber,
           Acetamiprid, 
           Clothianidin,
           Dinotefuran,
           Imidacloprid,
           Thiacloprid,
           Thiamethoxam) %>%
    gather(chemical, value, -site, -NWISRecordNumber) %>%
    mutate(site = zeroPad(site,8))
  
  just_neonic_remarks <- neonic %>%
    select(site = USGS.Site.ID,  NWISRecordNumber,
           R_Acetamiprid, 
           R_Clothianidin,
           R_Dinotefuran,
           R_Imidacloprid,
           R_Thiacloprid,
           R_Thiamethoxam
    ) %>%
    gather(chemical_rk, remark_cd, -site,  -NWISRecordNumber) %>%
    mutate(chemical = gsub("R_","",chemical_rk)) %>%
    select(-chemical_rk)%>%
    mutate(site = zeroPad(site,8))

  just_neonic <- left_join(just_neonic_data, 
                           just_neonic_remarks, by=c("site","chemical","NWISRecordNumber")) %>%
    left_join(select(tracking, NWISRecordNumber, pdate), by="NWISRecordNumber")

  just_NWIS <- select(NWIS, site=SiteID, pdate, pCode, value, remark_cd) %>%
    left_join(select(pCodeInfo, pCode=parameter_cd, chemical=casrn), by="pCode") 
  
  nwis_neonic <- bind_rows(just_neonic, just_NWIS)
  
  return(nwis_neonic)
}

remove_censor <- function(neonic_NWIS){
  
  neonic_NWIS$value[neonic_NWIS$remark_cd == "<"] <- 0
  neonic_NWIS$value[is.na(neonic_NWIS$value)] <- 0
  
  neonic_NWIS <- filter(neonic_NWIS, value != 0)
    
  return(neonic_NWIS)
  
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
                          Class = "Insecticide",stringsAsFactors = FALSE)
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
  
  chem_info <- select(chem_info, -class1, -pCode) %>%
    distinct()

  chem_info$`Chemical Name` <- gsub(", water, filtered, recoverable, nanograms per liter","",chem_info$`Chemical Name`)
  # Note, we're converting from ng to ug in the create_chemData function
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
  
create_toxExcel <- function(chem_data, chem_info, site_info, exclusions, file_out){
  

  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info,
                           "Exclude" = exclusions)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)

}

get_chem_sum <- function(chem_info, chem_data, site_info, exclusions){

  ACClong <- get_ACC(chem_info$CAS)
  ACClong <- remove_flags(ACClong)
  
  cleaned_ep <- clean_endPoint_info(endPointInfo)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(ACClong,
                                          filtered_ep,
                                          chem_data, 
                                          site_info, 
                                          chem_info,
                                          exclusions)
  return(chemicalSummary)
}

get_chem_bench <- function(benchmarks, chem_data, site_info, chem_info, exclusions){
  
  benchmarks <- benchmarks %>%
    rename(chnm = Compound,
           ACC_value = value) %>%
    filter(!is.na(CAS))
  
  filtered_ep <- select(benchmarks, endPoint) %>%
    distinct() %>%
    mutate(groupCol = "Aquatic Benchmark")
  
  chemicalSummary_bench <- get_chemical_summary(benchmarks,
                                                filtered_ep,
                                                chem_data, 
                                                site_info, 
                                                chem_info,
                                                exclusions)
  return(chemicalSummary_bench)
  
}

get_conc_summary <- function(chem_data, site_info, chem_info, exclusions){
  
  conc_ep <- select(chem_info, CAS, chnm=`Chemical Name`) %>%
    mutate(ACC_value = 1,
           endPoint = "Concentration") %>%
    filter(!is.na(CAS)) %>%
    distinct()
  
  filtered_ep <- select(conc_ep, endPoint) %>%
    distinct() %>%
    mutate(groupCol = "Concentrations")
  
  chemicalSummary_conc <- get_chemical_summary(conc_ep,
                                                filtered_ep,
                                                chem_data, 
                                                site_info, 
                                                chem_info,
                                                exclusions)
  return(chemicalSummary_conc)
  
}
