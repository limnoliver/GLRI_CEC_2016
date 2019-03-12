library(dplyr)
owc <- read.csv('owc_2016.csv', stringsAsFactors = FALSE, colClasses = c('pCode' = 'character'))

owc_merge <- select(owc, pCode) %>%
  unique() %>% 
  left_join(select(dataRetrieval::parameterCdFile, pCode = parameter_cd, parameter_nm)) %>%
  mutate(parameter_nm_short = gsub(', water, .+', '', parameter_nm)) %>%
  mutate(parameter_nm_short = gsub(', NWQL.+', '', parameter_nm_short))

owc_new <- owc %>%
  left_join(select(owc_merge, pCode, parameter_nm_short)) %>%
  mutate(srsname = ifelse(is.na(srsname), parameter_nm_short, srsname)) %>%
  select(-parameter_nm_short)

write.csv(owc_new, 'usgs_owc_2016.csv', row.names = FALSE)

pesticides <- read.csv('pesticides_2016.csv', stringsAsFactors = FALSE, colClasses = c('pCode' = 'character'))

pesticides_merge <- select(pesticides, pCode) %>%
  unique() %>% 
  left_join(select(dataRetrieval::parameterCdFile, pCode = parameter_cd, parameter_nm)) %>%
  mutate(parameter_nm_short = gsub(', water, .+', '', parameter_nm)) %>%
  mutate(parameter_nm_short = gsub(', NWQL.+', '', parameter_nm_short))

pesticides_new <- pesticides %>%
  left_join(select(pesticides_merge, pCode, parameter_nm_short)) %>%
  mutate(srsname = ifelse(is.na(srsname), parameter_nm_short, srsname)) %>%
  select(-parameter_nm_short)

write.csv(pesticides_new, 'usgs_pesticides_2016.csv', row.names = FALSE)
