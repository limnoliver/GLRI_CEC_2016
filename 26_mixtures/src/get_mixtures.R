# Mixtures
# Modifying this script from what Steve did for mixtures
# in general this script:

# 1. For each sample, sums EARs by endpoint
# 2. For each chemical, calculates percent contribution to each summed endpoint
# 3. Chooses a cutoff for percent contribution, and finds important mixtures.
# 4. Defines chemicals in important mixtures
sum_endpoints <- function(all_EARs, filter_cutoffs = TRUE) {
ear_cutoff <- 0.001

# note Steve limited each chemical to contributing to one endpoint by 
# only using the max EAR val for each site-date-chem combination
# I don't really know the rationale for this -- so am leaving out for now
# I've incorporated a slightly different filter below: I used only the max summed endpointEAR
# per site-date to represent sort of the "worst" mixture. Should touch base with Steve about this. 

if (filter_cutoffs) {
  summed_EARs <- all_EARs %>%
    group_by(site, shortName, date, endPoint) %>%
    mutate(sum_ear_endpoint = sum(EAR)) %>%
    ungroup() %>%
    mutate(chem_mix_contribution = (EAR/sum_ear_endpoint)*100) %>%
    filter(sum_ear_endpoint > ear_cutoff) %>%
    filter(chem_mix_contribution > 1)
} else {
  summed_EARs <- all_EARs %>%
    group_by(site, shortName, date, endPoint) %>%
    mutate(sum_ear_endpoint = sum(EAR)) %>%
    ungroup() %>%
    mutate(chem_mix_contribution = (EAR/sum_ear_endpoint)*100) %>%
    filter(chem_mix_contribution > 1)
}


return(summed_EARs)
}

calc_contr_chems <- function(summed_EARs) {
EAR_sum_endpoint <- summed_EARs %>%
  group_by(site, shortName, date, endPoint, sum_ear_endpoint) %>%
  summarize(n_contr_chems = n(),
            n_contr_parents = length(unique(parent_pesticide)),
            contr_chems = paste0(sort(chnm), collapse = ', '),
            contr_parents = paste0(sort(unique(parent_pesticide)), collapse = ', '),
            max_individual_contr = max(EAR)) %>%
  ungroup()
return(EAR_sum_endpoint)
}

calc_top_mixtures <- function(EAR_sum_endpoint, max_only = TRUE) {
  
  if (max_only) {
    top_mixtures <- EAR_sum_endpoint %>%
      group_by(site, shortName, date) %>%
      summarize(max_sum_ear_endpoint = max(sum_ear_endpoint),
                endPoint_top = endPoint[which.max(sum_ear_endpoint)],
                n_contr_chems = n_contr_chems[which.max(sum_ear_endpoint)],
                contr_chems = contr_chems[which.max(sum_ear_endpoint)][order(contr_chems[which.max(sum_ear_endpoint)])],
                max_individual_contr = max_individual_contr[which.max(sum_ear_endpoint)]) %>%
      mutate(prop_ind_contr = max_individual_contr/max_sum_ear_endpoint)
  } else {
    top_mixtures <- EAR_sum_endpoint %>%
      mutate(prop_ind_contr = max_individual_contr/sum_ear_endpoint,
             endPoint_top = endPoint,
             max_sum_ear_endpoint = sum_ear_endpoint)
  }
  
  return(top_mixtures)
}
# calculate metrics by chemical
# calculate the number of times a chemical was in a mixture
# exclude 1-compound mixtures first

calc_chem_mix_metrics <- function(top_mixtures, summed_EARs, out_file) {
  
top_mixes_2plus <- filter(top_mixtures, n_contr_chems >1)

top_mix_chems <- summed_EARs %>%
  left_join(select(top_mixes_2plus, site, shortName, date, endPoint = endPoint_top, prop_ind_contr)) %>%
  filter(!is.na(prop_ind_contr)) %>%
  group_by(site, shortName, date, endPoint, sum_ear_endpoint, prop_ind_contr, parent_pesticide) %>%
  summarize(chem_mix_contribution = sum(chem_mix_contribution)) %>%
  group_by(parent_pesticide) %>%
  summarize(times_in_mixes = n(),
            contribution_median = median(chem_mix_contribution),
            n_endpoints = length(unique(endPoint)),
            n_sites = length(unique(site))) %>%
  arrange(-times_in_mixes)

return(top_mix_chems)
}

summarize_mixtures <- function(top_mixtures) {
  mix_summary <- top_mixtures %>%
    filter(n_contr_chems >1) %>%
    group_by(contr_chems, n_contr_chems) %>%
    summarize(endPoint = paste(unique(endPoint_top), collapse = ', '), 
              mix_n_hits = n(),
              mix_n_hits_sites = length(unique(site)),
              mix_n_hits_samples = length(unique(paste0(site, date))),
              mix_n_hits_months = length(unique(lubridate::month(date))),
              mix_max_sum_EAR = max(max_sum_ear_endpoint),
              mix_median_sum_EAR = median(max_sum_ear_endpoint),
              mix_median_prop_individual_contr = median(prop_ind_contr))
    
}

summarize_by_n <- function(top_mixtures) {
  n_summary <-top_mixtures %>%
    filter(n_contr_chems >1) %>%
    group_by(n_contr_chems) %>%
    summarize(n_hits = n(),
              n_unique_mixes = length(unique(contr_chems)))
}

plot_mix_summary <- function(n_summary, mix_summary, top_mixtures, ear_sum, out_file) {
  p1 <- ggplot(n_summary, aes(x = n_contr_chems, y = n_hits)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label = n_unique_mixes), vjust = -0.25, color = 'red3') +
    scale_x_continuous(breaks = 2:8) +
    theme_bw() +
    labs(x = 'Number of chemicals in mixture', y = 'Hits (EARmix > 0.001)') +
    annotate('text', x = 6.6, y = 250, label = 'Number of \nunique mixtures', color = 'red3')
  
  top <- mix_summary %>% 
    filter(mix_n_hits_samples >= 10) %>%
    arrange(-mix_n_hits_samples)
    
  top_all <- filter(top_mixtures, contr_chems %in% top$contr_chems)  %>%
    mutate(top = TRUE)
  

  top_all_chems <- left_join(ear_sum, select(top_all, site, date, endPoint_top, contr_chems, top), by = c('site', 'date', 'endPoint' = 'endPoint_top')) %>%
    filter(!is.na(top))
  
  plot_summary <- top_all_chems %>%
    group_by(chnm, CAS, Class, contr_chems) %>%
    summarize(median_ear = median(EAR)) %>%
    left_join(top, by = 'contr_chems')
  
  plot_summary$contr_chems <- factor(plot_summary$contr_chems, levels = top$contr_chems)
  
  totals <- plot_summary %>%
    group_by(contr_chems, mix_n_hits_samples) %>%
    summarize(total = sum(median_ear)) %>%
    mutate(chnm = NA)
  
  plot_summary$chnm <- gsub('2,4-Dichlorophenoxyacetic acid', '2,4-D', plot_summary$chnm)
  p2 <- ggplot(plot_summary, aes(x = contr_chems, y = median_ear, fill = chnm)) +
    geom_bar(stat = 'identity', color = 'black') +
    geom_text(data = totals, aes(x = contr_chems, y = total, label = mix_n_hits_samples), vjust = -0.25, color = 'red3') +
    theme_bw() +
    geom_hline(yintercept = 0.001, linetype = 2) +
    theme(axis.text.x = element_blank(),
          legend.position = c(0.25, 0.7),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key.size = unit(0.3, "cm")) +
    labs(x = 'Unique mixtures', y = 'Median EARmix', fill = '') +
    annotate('text', x = 8, y = 0.016, label = 'Number of \nsample hits', color = 'red3')
  
  
  ggsave(out_file, cowplot::plot_grid(p1, p2), height = 3.6, width = 8)
    
  
  }

# sites - calculate max EARmix across samples, as well as the 
# number of months where there is an EARmix > 0.001 
calc_site_mix_metrics <- function(top_mixtures, all_samples) {
  all_dates <- all_samples %>%
    select(site, date) %>%
    distinct() %>%
    group_by(site) %>%
    summarize(n_samples = n())
    
  site_mix <- top_mixtures %>%
    left_join(all_dates) %>%
    group_by(site, shortName) %>%
    summarize(mix_hits_per_sample = round(n()/unique(n_samples), 1),
             mix_hit_n_endpoints = length(unique(endPoint_top)),
            mix_hit_n_months = length(unique(lubridate::month(date))),
            max_EARmix = round(max(max_sum_ear_endpoint), 4)) %>%
    arrange(-max_EARmix)
  
  if (length(unique(site_mix$site)) < length(unique(all_dates$site))) {
    missing_sites <- unique(all_dates$site)[!unique(all_dates$site) %in% unique(site_mix$site)]
    
    # these sites are missing because they had no hits
    max_value <- filter(EARsum_endpoint_nofilter, site %in% missing_sites) %>%
      group_by(site, shortName) %>%
      summarize(max_EARmix = round(max(sum_ear_endpoint), 4))
    
    dat_add <- data.frame(site = max_value$site, 
                          shortName = max_value$shortName,
                          mix_hits_per_sample = 0, 
                          mix_hit_n_endpoints = 0,
                          mix_hit_n_months = 0, 
                          max_EARmix = max_value$max_EARmix)
    
    site_mix <- bind_rows(site_mix, dat_add)
  }
    
    return(site_mix)
}
