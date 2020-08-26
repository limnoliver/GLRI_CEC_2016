# chems missing toxEval dat 
# for chems that are not in toxEval
# sites chems are present
# range of concentrations
# CAS num/chem name

tox_dat <- toxEval::tox_chemicals

# read in your own chems (has pcode but no other ID)
chem_dat <- make('pesticides')

# merge with pcode info (basically dataRetrieval::parameterCdFile)
pcode_info <- make('pCodeInfo')
classes <- make('classes')

chems <- left_join(chem_dat, pcode_info, by = c('pCode' = 'parameter_cd'))

# have to use longer parameter_nm because some shorter srsnames don't exist
missing_chems <- select(chems, parameter_nm,  casrn, pCode) %>%
  distinct() %>%
  filter(!casrn %in% unique(tox_dat$Substance_CASRN)) %>%
  filter(!is.na(pCode))

missing_chems_summary <- chems %>%
  filter(pCode %in% missing_chems$pCode) %>%
  group_by(pCode) %>%
  summarize(n_samples = n(),
            n_detect_samples = sum(!grepl('<', remark_cd)), 
            perc_sample_detect = round((n_detect_samples/n_samples)*100,0),
            n_sites = length(unique(SiteID)),
            n_detect_sites = length(unique(SiteID[!grepl('<', remark_cd)])),
            perc_site_detect = round((n_detect_sites/n_sites) *100,0),
            min_conc = ifelse(n_detect_samples == 0, NA, min(value[!grepl('<', remark_cd)])),
            max_conc = ifelse(n_detect_samples == 0, NA, max(value[!grepl('<', remark_cd)])),
            mean_conc = ifelse(n_detect_samples == 0, NA, mean(value[!grepl('<', remark_cd)])),
            median_conc = ifelse(n_detect_samples == 0, NA, median(value[!grepl('<', remark_cd)])))

missing_chems_all <- left_join(missing_chems, missing_chems_summary) %>%
  left_join(classes, by = c('casrn'='CAS'))

write.csv(missing_chems_all, 'pesticide_summary_no_toxeval.csv')
