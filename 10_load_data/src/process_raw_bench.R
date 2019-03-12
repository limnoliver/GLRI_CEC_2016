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