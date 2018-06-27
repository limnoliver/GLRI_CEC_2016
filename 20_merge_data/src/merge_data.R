library(dplyr)
library(tidyr)
library(readr)
library(dataRetrieval)
library(openxlsx)
library(toxEval)

merge_data <- function(tracking, pesticides_clean, neonics_clean, glyphosate_clean){

  tracking_sub <- filter(tracking, MediumCode %in% 'WS') %>%
    select(-MediumCode, -SampleTypeCode, -(ToxCast:Passive), -pdate)

  all_dat <- bind_rows(pesticides_clean, neonics_clean, glyphosate_clean) %>%
    filter(medium_cd %in% 'WS')
  
  all_dat <- left_join(all_dat, tracking_sub, by = c('SiteID', 'sample_dt' = 'Date'))
  
  return(all_dat)
}

remove_duplicate_chems <- function(merged_dat) {
  # a priori remove chemicals that were measured at same site-date
  # by two different methods -- e.g., imidacloprid, glyphosate
  # imidacloprid is duplicated in every instance where there was the pesticide schedule + neonics measured
  # per Michelle H.'s recommendation, we'll create a relationship between the two methods for imidacloprid
  # use the neonic value when available, and use the estimated value when just 2437 data is available
  
  imidacloprid <- filter(merged_dat, pCode %in% '68426') %>%
    select(SiteID, Site, source, pCode, sample_dt, value) %>%
    distinct() %>% # this is a temporary fix for the duplicated IHC sample
    spread(key = source, value = value)
  
  imidacloprid_rmk <- filter(merged_dat, pCode %in% '68426') %>%
    select(SiteID, Site, source, pCode, sample_dt, remark_cd) %>%
    distinct() %>% # this is a temporary fix for the duplicated IHC sample
    spread(key = source, value = remark_cd) %>%
    rename(neonic_rmk = neonic, pesticides_s2437_rmk = pesticides_s2437)
  
  imidacloprid <- left_join(imidacloprid, imidacloprid_rmk)
  
  imidacloprid_adl <- filter(imidacloprid, !is.na(neonic) & !is.na(pesticides_s2437)) %>%
    filter(!(neonic_rmk %in% '<') & !(pesticides_s2437_rmk %in% '<'))
  
  mod <- lm(imidacloprid_adl$neonic ~ imidacloprid_adl$pesticides_s2437)
  adjust_2437 <- as.numeric(mod$coefficients[2]) # on average, the neonic result is 0.48 that of s2437, 
  # so will adjust s2437 by multiplying by 0.48
  
  imidacloprid <- mutate(imidacloprid, combined = ifelse(!(is.na(neonic)), neonic, 
                                                         ifelse(pesticides_s2437_rmk %in% "<", pesticides_s2437, pesticides_s2437*adjust_2437)))
  
  imidacloprid <- mutate(imidacloprid, combined_rmk = ifelse(!is.na(neonic), neonic_rmk, pesticides_s2437_rmk))
  
  imidacloprid <- mutate(imidacloprid, source = ifelse(!is.na(neonic), 'neonic', 'pesticides_s2437'))
  imidacloprid <- select(imidacloprid, SiteID, sample_dt, source, pCode, combined, combined_rmk) %>%
    rename(value = combined, remark_cd = combined_rmk)
  
  # now need to remove all imidacloprid from merged_dat, and add "combined" value back in.
  fixed_imidacloprid <- left_join(imidacloprid, select(merged_dat, -value, -remark_cd), by = c('SiteID', 'sample_dt', 'pCode', 'source'))
  
  # now drop imidacloprid values that are above DL, and replace with fixed vals
  fixed_dat <- filter(merged_dat, !(pCode %in% '68426')) %>%
    bind_rows(fixed_imidacloprid)
    
  # drop the full glyphosate analysis + degradate, use the immunoassay
  # eventually should estimate degradate based on glyphosate ~ degradate relationship - see Mahler's work on topic
  glyphosate_drop <- c("62722", "62649")
  
  fixed_dat <- filter(fixed_dat, !(pCode %in% glyphosate_drop))
  
  # check fixed data for remaining duplicates
  
  # diagnostics:
  dup_chems <- filter(fixed_dat, !(SiteID %in% "04092750" & sample_dt %in% as.Date('2016-08-02'))) %>% # temp fix for IHC site on 8/2 duplicate data that is being fixed on NWIS
    group_by(SiteID, sample_dt, pCode) %>%
    filter(n()>1) 
  
  dup_pcodes <- filter(parameterCdFile, parameter_cd %in% unique(dup_chems$pCode))
  
  dup_chems_summary <- ungroup(dup_chems) %>%
    group_by(pCode) %>%
    select(pCode, SiteID, sample_dt) %>%
    distinct()
  
  # a single atrazine measurement randomly repeated on 6/7 from Maumee
  # two different values - nothing else repeated at that site/date
  # values are similar, just take first one, assume the second one is a replicate measure
  
  fixed_dat <- filter(fixed_dat, !(sample_dt %in% as.Date('2016-06-07') &
                                     sample_tm %in% '11:20' & 
                                     SiteID %in% '04193500' &
                                     pCode %in% '65065'))
  
  return(fixed_dat)

  
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
           Value = value,
           remark_cd,
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
           station_nm) %>%
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
