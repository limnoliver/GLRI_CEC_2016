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
  
  # also handle chemicals that have seperate pCodes but the same cas number. This will cause issues later as many
  # rely on unique casrns. 
  
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
  
  # now merge pcodes with cas numbers to find duplicate cas
  # check for duplicates
  dup_cas <- fixed_dat %>%
    left_join(select(parameterCdFile, pCode = parameter_cd, casrn, parameter_nm)) %>%
    filter(!is.na(casrn)) %>%
    select(pCode, casrn, parameter_nm) %>%
    distinct() %>%
    group_by(casrn) %>%
    summarize(count = n(), 
              pCode = paste(pCode, collapse = ', ')) %>%
    filter(count > 1)
  
  # only one casrn number is repeated
  # this is for Demethyl hexazinone B (68566)/Didemethyl hexazinone F (68574)
  # hex F has wrong cas number

  #fixed_dat <- filter(fixed_dat, pCode != '68574')
  # commenting out since we fix cas elsewhere
  
  return(fixed_dat)

  
}

calc_detect_limits <- function(reduced_dat){
  
  all_measured <- group_by(reduced_dat, pCode) %>%
    summarise(n_all = n())
  
  all_detected <- filter(reduced_dat, !(remark_cd %in% "<")) %>%
    group_by(pCode) %>%
    summarise(n_detect = n())
  
  all_bdl <- filter(reduced_dat, remark_cd %in% "<") %>%
    group_by(pCode) %>%
    summarise(n_bdl = n())
  
  bdls <- filter(reduced_dat, remark_cd %in% "<") %>%
    group_by(pCode) %>%
    summarize_at(vars(value), funs(mean, median, min, max))
  
  bdl_summ <- left_join(all_measured, all_detected) %>%
    left_join(all_bdl) %>%
    left_join(bdls) %>%
    mutate(n_detect = ifelse(is.na(n_detect), 0, n_detect)) %>%
    rename(value = median)
  
  bdl_summ$SiteID <- '01'
  bdl_summ$sample_dt <- '1986-11-13'
  bdl_summ$remark_cd <- NA
  
  
  return(bdl_summ)
}

remove_censor <- function(reduced_dat){
  
  reduced_dat$value[reduced_dat$remark_cd == "<"] <- 0
  #reduced_dat$value[is.na(reduced_dat$value)] <- 0
  
  # neonic_NWIS <- filter(neonic_NWIS, value != 0)
    
  return(reduced_dat)
  
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

create_chemData <- function(reduced_dat,  pCodeInfo){
  
  
  dat_with_cas <- left_join(reduced_dat, pCodeInfo, by = c('pCode' = 'parameter_cd')) 

    chem_data <- dat_with_cas %>%
    select(SiteID,
           `Sample Date` = sample_dt,
           Value = value,
           remark_cd,
           CAS = casrn,
           pCode, 
           parameter_units) %>%
    mutate(Value = if_else(parameter_units == "ng/l",Value/1000,Value)) %>%
    select(-parameter_units) %>%
    filter(!is.na(CAS)) # get rid of chemicals without CAS number - not sure if/where we should do this? 
  
  
  chem_data <- distinct(chem_data) # currently this just gets rid of duplicated IHC values on 8/2
  # change incorrect CAS numbers
  
  # Hexazinone TP F	68574 should have cas == 56611-55-3 from 56611-54-2
  chem_data$CAS[chem_data$pCode %in% '68574'] <- '56611-55-3'
  
  return(chem_data)
}

create_chemData2 <- function(reduced_dat,  pCodeInfo){
  
  dat_with_cas <- left_join(reduced_dat, pCodeInfo, by = c('pCode' = 'parameter_cd'))
  
  chem_data <- dat_with_cas %>%
    select(SiteID,
           `Sample Date` = sample_dt,
           Value = value,
           remark_cd,
           CAS = casrn,
           pCode, 
           parameter_units) %>%
    mutate(Value = if_else(parameter_units == "ng/l",Value/1000,Value)) %>%
    select(-parameter_units) # get rid of chemicals without CAS number - not sure if/where we should do this? 
  
  
  chem_data <- distinct(chem_data) # currently this just gets rid of duplicated IHC values on 8/2
  # change incorrect CAS numbers
  
  # Hexazinone TP F	68574 should have cas == 56611-55-3 from 56611-54-2
  chem_data$CAS[chem_data$pCode %in% '68574'] <- '56611-55-3'
  
  #chem_data
  
  
  
  return(chem_data)
}


find_missing_cas <- function(reduced_dat, chem_data) {
  chems_all <- unique(reduced_dat$pCode)
  chems_no_missing_cas <- unique(chem_data$pCode)
  
  chems_missing <- chems_all[-which(chems_all %in% chems_no_missing_cas)]
  return(chems_missing)
  

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
  
  chem_info <- select(chem_info, -class1) %>%
    distinct()

  chem_info$`Chemical Name` <- gsub(", water, filtered, recoverable, nanograms per liter","",chem_info$`Chemical Name`)
  
  # this all unnecessary now that we took care of duplicate methods/measures
  # chem_info$`Chemical Name`[chem_info$CAS == "138261-41-3_GLRI"] <- "Imidacloprid"
  # chem_info$Class[chem_info$CAS == "138261-41-3_GLRI"] <- "Insecticide"
  # 
  # chem_info$`Chemical Name`[chem_info$CAS == "56611-54-2_68574"] <- "Didemethyl hexazinone F"
  # chem_info$Class[chem_info$CAS == "56611-54-2_68574"] <- "Deg - Herbicide"
  # 
  chem_info$`Chemical Name`[chem_info$CAS == "1066-51-9"] <- "Aminomethylphosphonic acid"
  chem_info$`Chemical Name`[chem_info$CAS == '2303-17-5'] <- 'Triallate'
  chem_info$`Chemical Name`[chem_info$CAS == '94-75-7'] <- '2,4-D'
  
  
  # this needs to happen because glyphosate is in micrograms per liter. 
  # be sure this gets taken care of somewhere. 
  
  chem_info$`Chemical Name`[chem_info$CAS == "1071-83-6"] <- "Glyphosate"
  # chem_info$`Chemical Name`[chem_info$CAS == "1071-83-6_99960"] <- "Glyphosate"
  
  chem_info <- distinct(chem_info)
  
  # Note, we're converting from ng to ug in the create_chemData function
  return(chem_info)
}

fix_missing_cas <- function(chem_info_complete, chem_crosswalk) {
  
  info <- left_join(chem_info_complete, chem_crosswalk) %>%
    mutate(parent_pesticide = ifelse(is.na(parent_pesticide), `Chemical Name`, parent_pesticide))
  
  missing <- filter(info, is.na(Class))
  
  info$Class[info$pCode %in% '68561'] <- 'Deg - Insecticide'
  info$Class[info$pCode %in% '68685'] <- 'Deg - Herbicide'
  info$Class[info$pCode %in% c('68623', '68675')] <- 'Deg - Insecticide'
  info$Class[info$pCode %in% c('68581', '68583')] <- 'Deg - Herbicide'
  info$Class[info$pCode %in% '68563'] <- 'Deg - Insecticide'
  info$Class[info$pCode %in% c('68617', '68619', '68620')] <- 'Deg - Herbicide'
  info$Class[info$pCode %in% '68713'] <- 'Deg - Herbicide'
  info$Class[info$pCode %in% '68690'] <- 'Deg - Herbicide'
  info$Class[info$pCode %in% '68694'] <- 'Deg - Insecticide'
  info$Class[info$pCode %in% c('68575', '68714')] <- 'Deg - Herbicide'
  
  # additionally, there are issues with some mismatching compound to parent names
  info$parent_pesticide[info$parent_pesticide %in% 'Bentazone'] <- 'Bentazon'
  
  # add parent compounds that weren't measured but are used for degradate matching
  info <- add_row(info, 
                  `Chemical Name` = c('Chlorothalonil', 'Dacthal', 'Benomyl', 'Parathion-ethyl', 'Parathion-methyl'), 
                  Class = c('Fungicide', 'Herbicide', 'Fungicide', 'Insecticide', 'Insecticide'),
                  parent_pesticide = c('Chlorothalonil', 'Dacthal', 'Benomyl', 'Parathion-ethyl', 'Parathion-methyl'), 
                  compound = c('Chlorothalonil', 'Dacthal', 'Benomyl', 'Parathion-ethyl', 'Parathion-methyl'), 
                  CAS = c('1897-45-6', '1861-32-1', '17804-35-2', '56-38-2', '298-00-0'), 
                  MlWt = c(265.90, 331.95, 290.32, 291.26, 263.2),
                  pCode = paste0('fake_', 2:6))
  
  # fix glyphosate and deg
  info$Class[info$pCode %in% '62649'] <- 'Deg - Herbicide'
  info$parent_pesticide[info$pCode %in% '62649'] <- 'Glyphosate'
  info$MlWt[info$pCode %in% '62649'] <- 111.037
  info$MlWt[info$pCode %in% '99960'] <- 169.073
  
  
  # complete permethrin data
  info$Class[info$CAS %in% '52645-53-1'] <- 'Insecticide'
  info$compound[info$CAS %in% '52645-53-1'] <- 'Permethrin'
  info$MlWt[info$CAS %in% '52645-53-1'] <- 391.29
  
  # add CAS numbers to those compounds missing CAS numbers
  missing_cas <- filter(info, is.na(CAS)) %>%
    mutate(CAS = paste0('fakecas_', 1:n()))
  
  info <- filter(info, !is.na(CAS)) %>%
    bind_rows(missing_cas)
  
  return(info)
  
}

complete_cas <- function(chem_dat_complete, chem_master, filter_dl) {

  
  out <- chem_dat_complete %>%
    left_join(select(chem_master, CAS2 = CAS, pCode)) %>%
    mutate(CAS = ifelse(is.na(CAS), CAS2, CAS)) %>%
    select(-CAS2)
  
  if(filter_dl) {
    out <- filter(out, !remark_cd %in% '<')
    
  }
  
  
  return(out)
  
}
create_tox_siteInfo <- function(sites){

  siteInfo <- sites %>%
    select(SiteID = site_no,
           `Short Name` = shortName,
           dec_lat = dec_lat_va,
           dec_lon = dec_long_va,
           map_nm,
           station_nm) %>%
    mutate(`Short Name` = as.character(`Short Name`)) %>%
    mutate(site_grouping = "All")
  
  siteInfo$`Short Name`[siteInfo$SiteID %in% '04027000'] <- 'Bad'
  siteInfo$`Short Name`[siteInfo$SiteID %in% '04199500'] <- 'Vermillion'
  siteInfo$`Short Name`[siteInfo$SiteID %in% '04249000'] <- 'Oswego'
  
  siteInfo <- select(siteInfo, -map_nm, -station_nm)
  
  
  siteInfo$site_grouping[siteInfo$`Short Name` %in%  c("Bad","StLouis")] <- "Lake Superior"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Fox","Manitowoc","Milwaukee",
                                                      "IndianaHC","StJoseph","GrandMI")] <- "Lake Michigan"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Saginaw")] <- "Lake Huron"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Rouge","Clinton","Maumee","Vermillion","Cuyahoga")] <- "Lake Erie"
  siteInfo$site_grouping[siteInfo$`Short Name` %in% c("Genesee","Oswego")] <- "Lake Ontario"
  
  return(siteInfo)
  
}
  
create_toxExcel <- function(chem_data, chem_info, site_info, exclusions, file_out){
  
  
  list_of_datasets <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = site_info,
                           "Exclude" = exclusions)
  write.xlsx(list_of_datasets, file = file_out, append=TRUE)

}

create_WQExcel <- function(chem_data, chem_info, site_info, exclusions, benchmarks, file_out){
 
  benchmarks_new <- data.frame(CAS = chem_info$CAS, 
                           orig_name = chem_info$`Chemical Name`,
                           stringsAsFactors = FALSE) %>%
    left_join(select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by="CAS") %>%
    left_join(select(benchmarks, CAS, endPoint, Value = value), by="CAS") %>%
    mutate(groupCol = "WQ") %>%
    filter(!is.na(Value))

  nrow(benchmarks)
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
  #chem_data$CAS[chem_data$CAS == "138261-41-3_GLRI"] <- "138261-41-3"
  
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
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(tox_list$chem_info$CAS)

  ACClong <- remove_flags(ACClong)
  
  cleaned_ep <- clean_endPoint_info(toxEval::end_point_info)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
  return(chemicalSummary)
}

get_chem_bench <- function(data_file){
  
  tox_list <- create_toxEval(data_file)
  
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  ncol(tox_list)

  chemicalSummary_bench <- get_chemical_summary(tox_list)
  
  return(chemicalSummary_bench)
  
}

get_conc_summary <- function(data_file){
  
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)

  chemicalSummary_conc <- get_chemical_summary(tox_list)
  return(chemicalSummary_conc)
  
}
