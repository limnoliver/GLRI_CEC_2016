# Mixtures
# Modifying this script from what Steve did for mixtures
# in general this script:

# 1. For each sample, sums EARs by endpoint
# 2. For each chemical, calculates percent contribution to each summed endpoint
# 3. Chooses a cutoff for percent contribution, and finds important mixtures.
# 4. Defines chemicals in important mixtures

all_EARs <- make('chemicalSummary')
ear_cutoff <- 0.001

# note Steve limited each chemical to contributing to one endpoint by 
# only using the max EAR val for each site-date-chem combination
# I don't really know the rationale for this -- so am leaving out for now
# I've incorporated a slightly different filter below: I used only the max summed endpointEAR
# per site-date to represent sort of the "worst" mixture. Should touch base with Steve about this. 
summed_EARs <- all_EARs %>%
  group_by(site, shortName, date, endPoint) %>%
  mutate(sum_ear_endpoint = sum(EAR)) %>%
  ungroup() %>%
  mutate(chem_mix_contribution = (EAR/sum_ear_endpoint)*100) %>%
  filter(sum_ear_endpoint > ear_cutoff) %>%
  filter(chem_mix_contribution > 1) 

EAR_sum_endpoint <- summed_EARs %>%
  group_by(site, shortName, date, endPoint, sum_ear_endpoint) %>%
  summarize(n_contr_chems = n(),
            contr_chems = paste0(chnm, collapse = ', '),
            max_individual_contr = max(EAR)) %>%
  ungroup()

top_mixtures <- EAR_sum_endpoint %>%
  group_by(site, shortName, date) %>%
  summarize(max_sum_ear_endpoint = max(sum_ear_endpoint),
            endPoint_top = endPoint[which.max(sum_ear_endpoint)],
            n_contr_chems = n_contr_chems[which.max(sum_ear_endpoint)],
            contr_chems = contr_chems[which.max(sum_ear_endpoint)][order(contr_chems[which.max(sum_ear_endpoint)])],
            max_individual_contr = max_individual_contr[which.max(sum_ear_endpoint)]) %>%
  mutate(prop_ind_contr = max_individual_contr/max_sum_ear_endpoint)

# calculate metrics by chemical
# calculate the number of times a chemical was in a mixture
# exclude 1-compound mixtures first

top_mixes_2plus <- filter(top_mixtures, n_contr_chems >1)

top_mix_chems <- summed_EARs %>%
  left_join(select(top_mixes_2plus, site, shortName, date, endPoint = endPoint_top, prop_ind_contr)) %>%
  filter(!is.na(prop_ind_contr)) %>%
  group_by(chnm) %>%
  summarize(times_in_mixes = n(),
            contribution_median = median(chem_mix_contribution))

# sites - calculate max EARmix across samples, as well as the 
# number of months where there is an EARmix > 0.001 
site_mix <- top_mixtures %>%
  group_by(site, shortName) %>%
  summarize(n_mix_hit_months = length(unique(lubridate::month(date))),
            max_EARmix = max(max_sum_ear_endpoint))

