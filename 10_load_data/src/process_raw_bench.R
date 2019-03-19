process_benchmarks <- function(in_file, out_file) {

 # raw benchmark data downloaded from https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/aquatic-life-benchmarks-and-ecological-risk
 bench <- read.csv(in_file, stringsAsFactors = FALSE, 
                   na.strings = c('', 'NA'), strip.white = TRUE, 
                   col.names = c('Compound', 'Year', 'CAS', 'AcuteFish', 'ChronicFish',
                                 'AcuteInvert', 'ChronicInvert', 'AcuteNonvascularPlants', 
                                 'AcuteVascularPlants', 'MaxConcAquaticLifeCriteriaCMC',  'ContinuousConcAquaticLifeCriteriaCMC'))
 bench_long <- bench %>% 
   gather(key = 'endPoint', value = 'value', AcuteFish:ContinuousConcAquaticLifeCriteriaCMC) %>%
   mutate(remark = case_when(
     grepl('<', value) ~ '<',
     grepl('>', value) ~ '>')) %>%
   mutate(value = gsub('<|> ', '', value)) %>%
   select(Compound, Year, CAS, endPoint, value, remark)
 
 write.csv(bench_long, file = out_file, row.names = FALSE)
 
}

process_hsbls <- function(in_pesticides, in_degradates, out_file) {
   pest <- read.csv(in_pesticides, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
   deg <- read.csv(in_degradates, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
   
   all <- bind_rows(pest, deg)
   
   names(all) <- c('Compound', 'CAS', 'pCode', 'ChemClass', 'EPAMaximumContaminantLevels', 'EPAChronicNoncancerHHBP',
                   'EPACarcinogenicHHBP', 'USGSChronicNoncancerHBSL', 'USGSCancerHBSL', 'remark')
   
   all_long <- all %>%
      tidyr::gather(key = 'endPoint', value = 'value', EPAMaximumContaminantLevels:USGSCancerHBSL) %>%
      filter(!is.na(value)) %>%
      filter(grepl('Noncancer', endPoint)) %>%
      mutate(value = as.numeric(value)) %>%
      select(Compound, CAS, endPoint, value, remark)
   
   write.csv(all_long, file = out_file, row.names = FALSE)
}