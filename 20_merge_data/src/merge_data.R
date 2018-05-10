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
  
  # neonic_NWIS$value[neonic_NWIS$remark_cd == "<"] <- 0
  # neonic_NWIS$value[is.na(neonic_NWIS$value)] <- 0
  
  # neonic_NWIS <- filter(neonic_NWIS, value != 0)
    
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
    mutate(Value = if_else(is.na(units) | units == "ng/l",Value/1000,Value) ) %>%
    select(-units)

  chem_data$CAS[!is.na(chem_data$pCode) & chem_data$pCode == "68574"] <- "56611-54-2_68574"
  # chem_data$CAS[!is.na(chem_data$pCode) & chem_data$pCode == "99960"] <- "1071-83-6_99960"
  
  chem_data$CAS[chem_data$CAS == "Acetamiprid"] <- "135410-20-7"
  chem_data$CAS[chem_data$CAS == "Thiamethoxam"] <- "153719-23-4"
  chem_data$CAS[chem_data$CAS == "Thiacloprid"] <- "111988-49-9"
  chem_data$CAS[chem_data$CAS == "Imidacloprid"] <- "138261-41-3_GLRI"
  chem_data$CAS[chem_data$CAS == "Dinotefuran"] <- "165252-70-0"
  chem_data$CAS[chem_data$CAS == "Clothianidin"] <- "210880-92-5"
  
  chem_data <- select(chem_data, -casrn, -Class) %>%
    distinct()

  return(chem_data)
  
  
}


create_tox_chemInfo <- function(chem_data, special_cas, pCodeInfo, classes){


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
  
  chem_info$`Chemical Name`[chem_info$CAS == "138261-41-3_GLRI"] <- "Imidacloprid"
  chem_info$Class[chem_info$CAS == "138261-41-3_GLRI"] <- "Insecticide"
  
  chem_info$`Chemical Name`[chem_info$CAS == "56611-54-2_68574"] <- "Didemethyl hexazinone F"
  chem_info$Class[chem_info$CAS == "56611-54-2_68574"] <- "Deg - Herbicide"
  
  chem_info$`Chemical Name`[chem_info$CAS == "1066-51-9"] <- "Aminomethylphosphonic acid"
  chem_info$`Chemical Name`[chem_info$CAS == "1071-83-6"] <- "Glyphosate"
  # chem_info$`Chemical Name`[chem_info$CAS == "1071-83-6_99960"] <- "Glyphosate"
  
  chem_info <- distinct(chem_info)
  
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
  
  siteInfo$site_grouping[siteInfo$`Short Name` %in%  c("Bad","StLouis")] <- "Lake Superior"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Fox","Manitowoc","Milwaukee",
                                                      "IndianaHC","StJoseph","GrandMI")] <- "Lake Michigan"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Saginaw")] <- "Lake Huron"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Rouge","Clinton","Maumee","Vermillion","Cuyahoga")] <- "Lake Erie"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Genesee","Oswego")] <- "Lake Ontario"
  
  return(siteInfo)
  
}
  
create_toxExcel <- function(chem_data, chem_info, site_info, exclusions, file_out){
  
  chem_data$CAS[chem_data$CAS == "56611-54-2_68574"] <- "56611-54-2"
  # chem_data$CAS[chem_data$CAS == "1071-83-6_99960"] <- "1071-83-6"
  chem_data$CAS[chem_data$CAS == "138261-41-3_GLRI"] <- "138261-41-3"

  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info,
                           "Exclude" = exclusions)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)

}

create_WQExcel <- function(chem_data, chem_info, site_info, exclusions, benchmarks, file_out){
  
  chem_data$CAS[chem_data$CAS == "56611-54-2_68574"] <- "56611-54-2"
  # chem_data$CAS[chem_data$CAS == "1071-83-6_99960"] <- "1071-83-6"
  chem_data$CAS[chem_data$CAS == "138261-41-3_GLRI"] <- "138261-41-3"
  
  benchmarks_new <- data.frame(CAS = chem_info$CAS, 
                           orig_name = chem_info$`Chemical Name`,
                           stringsAsFactors = FALSE) %>%
    left_join(select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by="CAS") %>%
    left_join(select(benchmarks, CAS, endPoint, Value = value), by="CAS") %>%
    mutate(groupCol = "WQ")
  
  benchmarks_new$chnm[is.na(benchmarks_new$chnm)] <- benchmarks_new$orig_name[is.na(benchmarks_new$chnm)]
  
  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info,
                           "Exclude" = exclusions,
                           "Benchmarks" = benchmarks_new)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)
  
}

create_ConcExcel <- function(chem_data, chem_info, site_info, exclusions, file_out){
  
  chem_data$CAS[chem_data$CAS == "56611-54-2_68574"] <- "56611-54-2"
  # chem_data$CAS[chem_data$CAS == "1071-83-6_99960"] <- "1071-83-6"
  chem_data$CAS[chem_data$CAS == "138261-41-3_GLRI"] <- "138261-41-3"
  
  benchmarks <- data.frame(CAS = chem_info$CAS, 
                           orig_name = chem_info$`Chemical Name`,
                           stringsAsFactors = FALSE) %>%
    left_join(select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by="CAS") %>%
    mutate(endPoint = "Concentration",
           Value = 1,
           groupCol = "Concentration")
  
  benchmarks$chnm[is.na(benchmarks$chnm)] <- benchmarks$orig_name[is.na(benchmarks$chnm)]

  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info,
                           "Exclude" = exclusions,
                           "Benchmarks" = benchmarks)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)
  
}

get_chem_sum <- function(data_file){

  tox_list <- create_toxEval(data_file)
  tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(tox_list$chem_info$CAS)
  ACClong <- remove_flags(ACClong)
  
  cleaned_ep <- clean_endPoint_info(endPointInfo)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
  return(chemicalSummary)
}

get_chem_bench <- function(data_file){
  
  tox_list <- create_toxEval(data_file)
  tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)

  chemicalSummary_bench <- get_chemical_summary(tox_list)
  
  return(chemicalSummary_bench)
  
}

get_conc_summary <- function(data_file){
  
  tox_list <- create_toxEval(data_file)
  tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)

  chemicalSummary_conc <- get_chemical_summary(tox_list)
  return(chemicalSummary_conc)
  
}
