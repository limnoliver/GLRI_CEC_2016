
#Load toxEval files and generate chemical summary for surface water pesticides

# This is where the data are located on Luke's USGS computer
# path_to_data <- c('C:/Users/lloken/DOI/Corsi, Steven R - GLRI CECs/2016/Manuscripts/Pesticides_passive')


#Choose date filer
# date_filter <-as.Date(range(c(sites_2016$Date_in, sites_2016$Date_out)))
date_filter <- as.Date(c("2016-06-01", "2016-07-21")) #Passive deploy dates

tox_list_surface <- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides.xlsx")))
# tox_list_surface$chem_site$site_grouping <- factor(tox_list_surface$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))

#passive tox_list
tox_list_passive<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/Passive2016_ToxEval.xlsx")))
tox_list_passive$chem_site$site_grouping <- factor(tox_list_passive$chem_site$site_grouping, c('MN', 'WI', 'IL', 'IN', 'MI', 'OH', 'NY'))


#Combine exclusions
exclusions <- tox_list_passive$exclusions %>%
  # dplyr::select(-chnm) %>%
  bind_rows(tox_list_surface$exclusions) %>%
  distinct()

tox_list_surface$exclusions <- exclusions

#Rename to match passive
tox_list_surface$chem_site <- tox_list_surface$chem_site %>%
  rename(Lake = site_grouping) %>%
  mutate(Lake = gsub("Lake ", "", Lake)) %>%
  mutate(`Short Name` = gsub('St', 'St. ', `Short Name`)) %>%
  mutate(`Short Name` = gsub('GrandMI', 'Grand', `Short Name`)) %>%
  mutate(`Short Name` = gsub('IndianaHC', 'Indiana Harbor', `Short Name`)) %>%
  mutate(`Short Name` = gsub('Vermillion', 'Vermilion', `Short Name`)) %>%
  left_join(tox_list_passive$chem_site[,c("SiteID", 'site_grouping')]) 

#Manually change chemical classes and State IDs
tox_list_surface$chem_info$Class[tox_list_surface$chem_info$CAS == '78-48-8'] = 'Herbicide'

tox_list_surface$chem_site$site_grouping[which(tox_list_surface$chem_site$SiteID == "04157000")] <- 'MI'

tox_list_surface$chem_site$site_grouping[which(tox_list_surface$chem_site$SiteID == "04085427")] <- 'WI'

#Flag some ACCs
ACClong <- get_ACC(tox_list_surface$chem_info$CAS)
# ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical","ACCLessThan"))
ACClong <- remove_flags(ACClong)
cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             remove_groups = c('Background Measurement','Undefined'))

chemicalSummary_surface <- get_chemical_summary(tox_list_surface, 
                                        ACClong, 
                                        filtered_ep)
#Add groups to summary
chemicalSummary_surface <- tox_list_surface$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>%
  right_join(chemicalSummary_surface) 

#Change saginaw river site to match pocis
chemicalSummary_surface$site[which(chemicalSummary_surface$site=='04157000')] <- "04157005"


#PLot order
surface_ID_order <- tox_list_surface$chem_site$SiteID

surface_ID_order[which(surface_ID_order=='04157000')] <- "04157005"


##################################
# custom benchmarks
# same file as select data, but use custom benchmarks
# use this for toxicity quotient (TQ)

tox_bench_list_surface<- create_toxEval(file_in(file.path(path_to_data, "ToxEvalFiles/WQ_pesticides_Bench.xlsx")))

#replace chem data on benchmarks file with the EAR file
#Note this also appplies any filters added above
tox_bench_list_surface$chem_data <- tox_list_surface$chem_data
tox_bench_list_surface$chem_info <- tox_list_surface$chem_info

chemicalSummary_bench_surface <- get_chemical_summary(tox_bench_list_surface)

chemicalSummary_bench_surface <- tox_bench_list_surface$chem_site %>%
  rename(site = SiteID,
         shortName = `Short Name`) %>% 
  right_join(chemicalSummary_bench_surface) 


#Change saginaw river site to match pocis
chemicalSummary_bench_surface$site[which(chemicalSummary_bench_surface$site=='04157000')] <- "04157005"



