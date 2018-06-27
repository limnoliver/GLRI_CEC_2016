# concentration summaries by site
# sum of compounds
library(dplyr)
library(ggplot2)
library(lubridate)

#conc_dat <- make('chemicalSummary_conc')
all_chems <- as.character(unique(conc_dat$chnm))
all_casrn <- as.character(unique(conc_dat$CAS))
tracking <- make('tracking')
pesticides <- make('pesticides')
neonic <- make('neonic')
glyphosate <- make('glyphosate')

merged_dat <- make('merged_dat')

# calculate how many compounds we should expect from NWIS pull
schedule_pCodes <- make('schedule_pCodes') # 249 compounds in schedule
exclude_pCodes <- make('pCodesExclude') # looks like all surrogates removed
nwis_missing <- schedule_pCodes[which(!(schedule_pCodes$`Parameter Code` %in% unique(NWIS$pCode))), ] # 24 compounds missing from NWIS compared to schedule - okay, these are all surrogates and "Sample Volume"
length(unique(NWIS$pCode)) #226 compounds, which is correct because 249 - 24 missing compounds = 225 compounds + glyphosate which added to pcode list
# so when just pesticides = 225, glyphosate can be 1 or 3, neonics 6
# possible options 225, 226, 228, 231, 232, 234

test <- group_by(merged_dat, SiteID, sample_dt) %>%
  summarize(n = n())

nwis_pcodes <- nrow(schedule_pCodes) - length(exclude_pCodes)

## neonic summary
neonic_sum <- neonic %>%
  mutate(date = date(pdate)) %>%
  group_by(SiteID, date) %>%
  summarize(n = n())

neonic_sum

## river rouge missing neonic data
rr_neonic <- neonic %>%
  mutate(date = date(pdate)) %>%
  filter(USGS.Site.ID == '04166500' & date == as.Date('2016-06-14'))

# number of chemicals just from NWQL schedule
# glyphosate = 65065
all_nwis_pcodes <- unique(NWIS$pCode)
n_NWIS <- NWIS %>%
  group_by(SiteID, pdate) %>%
  summarize(n_chem = n(), 
            n_unique_chem = length(unique(pCode)))

one_chem <- filter(n_NWIS, n_chem == 1)
one_chem_dat <- left_join(one_chem, NWIS)
one_chem_pcodes <- unique(one_chem_dat$pCode)
one_chem_names <- dataRetrieval::parameterCdFile %>%
  filter(parameter_cd %in% one_chem_pcodes) %>%
  select(parameter_nm)
n_NWIS
# check if output matches tracking to identify gaps

test <- conc_dat %>%
  mutate(date_no_time = date(date)) %>%
  group_by(site, date_no_time) %>%
  summarize(neonic = ifelse('Clothianidin' %in% chnm, TRUE, FALSE),
            n_glyphosate = length(grep('Glyphosate', chnm)), 
            n_chems = n(),
            n_dup_chems = n()-length(unique(chnm)), 
            n_dup_casrn = n() - length(unique(CAS)),
            n_missing_chems = length(which(!(all_chems %in% chnm))),
            names_missing_chems = paste(all_chems[which(!(all_chems %in% chnm))], collapse = ', ')) %>%
  arrange(site, date_no_time)

dup_chems <- conc_dat %>%
  mutate(date_no_time = date(date)) %>%
  group_by(site, date_no_time) %>%
  filter(duplicated(chnm)) %>%
  summarize(dup_chems = paste(chnm, collapse = ', '))

test_tracking <- tracking %>%
  select(Site, SiteID, Date, Time, NWISRecordNumber, Glyphosate, Neonics, ToxCast, Comments, pdate_new, pdate) %>%
  left_join(test, by = c('SiteID' = 'site', 'Date' = 'date_no_time')) %>%
  left_join(dup_chems, by = c('SiteID' = 'site', 'Date' = 'date_no_time')) %>%
  mutate(neonic_problem = ifelse((is.na(Neonics) & neonic == FALSE)|(!is.na(Neonics) & neonic == TRUE), FALSE, TRUE),
         glyphosate_problem = ifelse((is.na(Glyphosate) & n_glyphosate %in% 0)|(!is.na(Glyphosate) & n_glyphosate %in% 1), FALSE, TRUE))

# Demethyl hexazinone B (56611-54-2) is an oft repeated chemical - figure out what is going on here
dhb <- filter(conc_dat, chnm == 'Demethyl hexazinone B')
dhb_nwis <- filter(NWIS, pCode == '68566')
dhf_nwis <- filter(NWIS, pCode == '68574')
head(NWIS)
# 215 unique chemicals - probably 217 because glysophate is measured
# two different ways, and chemical overlap between neonics and schedule

# look at glyphosate - measured two ways
glyphosate <- conc_dat %>%
  filter(chnm == "Glyphosate")


# find which 


 
 ##
 # glyphosate = 62722
 # 99960 = immuno
 
 pcodes <- parameterCdFile
View(pcodes[grep('glyphosate', pcodes$parameter_nm, ignore.case = T), ])

# Demethyl hexazinone
dh <- grep('hexazinone', ACC$chnm, ignore.case = T)
View(ACC[dh, ]) # not in toxEval
benchmarks[grep('hexazinone', benchmarks$Compound, ignore.case = T),]
pcodes[grep('hexazinone', pcodes$parameter_nm), ]
