# Mixtures
# Modifying this script from what Steve did for mixtures
# in general this script:

# 1. For each sample, sums EARs by endpoint
# 2. For each chemical, calculates percent contribution to each summed endpoint
# 3. Chooses a cutoff for percent contribution, and finds important mixtures.
# 4. Defines chemicals in important mixtures

sum_by_endpoint <- function(all_EARs) {
  
  summed_EARs <- all_EARs %>%
    group_by(site, shortName, date, endPoint) %>%
    summarize(sum_ear_endpoint = sum(EAR),
              max_ear_contr = max(EAR))
  
  site_max_individual <- all_EARs %>%
    group_by(site, date) %>%
    summarize(length(which(EAR>0.001)))
  
  site_hits <- group_by(summed_EARs, site, date) %>%
    summarize(max_endpoint_sum = max(sum_ear_endpoint)) %>%
    group_by(site) %>%
    summarize(n_hits = length(which(max_endpoint_sum > 0.001)))
  
  # of the hits, how many would not have been hits without the mixtures?
  mix_required <- summed_EARs %>%
    filter(sum_ear_endpoint > 0.001) %>%
    filter(max_ear_contr <= 0.001) %>%
    mutate(mix_required = TRUE)
  
  # how many individual compound hits are excluded from the 1% rule?
  individual <- all_EARs %>%
    group_by(site, shortName, date, endPoint) %>%
    mutate(sum_ear_endpoint = sum(EAR),
           individual_contr = (EAR/sum_ear_endpoint)*100) %>%
    filter(EAR > 0.001 & individual_contr < 1)
  
  # how much does the EARmix increase when considering mixtures?

  max_mix <- summed_EARs %>%
    mutate(increase = ((sum_ear_endpoint - max_ear_contr)/max_ear_contr)*100) %>%
    filter(increase > 0)
  
  site_increase <- group_by(max_mix, site) %>%
    summarize(median = median(increase),
              mean = mean(increase),
              n_increase = sum(increase > 1),
              n = n())
  
  detects <- chem_conc %>%
    group_by(site) %>%
    summarize(n_detect = length(unique(CAS)))
  
  site_increase <- left_join(site_increase, detects)
  
  # test non degradate
  hits_nodegs <- chem_ear %>%
    group_by(site, date) %>%
    summarize(hit_endpoint = length(which(EAR > 0.001)))
  
  hits_earmix <- chem_ear %>%
    group_by(site, date, endPoint) %>%
    summarize(sumEAR = sum(EAR)) %>%
    group_by(site, date) %>%
    summarize(hits = length(which(sumEAR > 0.001)))
  
  
  site_mix <- summed_EARs %>%
    filter(sum_ear_endpoint > 0.001) %>%
    left_join(select(mix_required, site, date, endPoint, mix_required)) %>%
    mutate(mix_required = ifelse(is.na(mix_required), FALSE, mix_required)) %>%
    group_by(site, shortName, date) %>%
    summarize(hits = n(), mix_dependent = sum(mix_required))
  
  head(chem_ear)

  
  
  return(summed_EARs)
  
}

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
    filter(!(chem_mix_contribution < 1 & EAR < 0.001))
} else {
  summed_EARs <- all_EARs %>%
    group_by(site, shortName, date, endPoint) %>%
    mutate(sum_ear_endpoint = sum(EAR)) %>%
    ungroup() %>%
    mutate(chem_mix_contribution = (EAR/sum_ear_endpoint)*100) %>%
    filter(!(chem_mix_contribution < 1 & EAR < 0.001))
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
    #filter(n_contr_chems >1) %>%
    group_by(contr_chems, n_contr_chems, n_contr_parents, contr_parents) %>%
    summarize(endPoint = paste(unique(endPoint_top), collapse = ', '), 
              mix_n_hits_endpoints = length(unique(endPoint_top)),
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
    #filter(n_contr_chems >1) %>%
    group_by(n_contr_chems) %>%
    summarize(n_hits = n(),
              n_unique_mixes = length(unique(contr_chems)))
}

plot_mix_summary <- function(n_summary, mix_summary, top_mixtures, ear_sum, chem_master, out_file) {
  p1 <- ggplot(n_summary, aes(x = n_contr_chems, y = n_hits)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label = n_unique_mixes), vjust = -0.25, color = 'red3', size = 3) +
    scale_x_continuous(breaks = 1:15) +
    theme_bw() +
    labs(x = 'Number of chemicals in mixture', y = 'Hits (EARmix > 0.001)') +
    annotate('text', x = 6.7, y = 550, label = 'Number of \nunique mixtures', color = 'red3', size = 3)
  
  top <- mix_summary %>% 
    filter(n_contr_chems > 1) %>%
    filter(mix_n_hits_samples > 10) %>%
    arrange(-mix_n_hits_samples)
    
  top_all <- filter(top_mixtures, contr_chems %in% top$contr_chems)  %>%
    mutate(top = TRUE)
  

  top_all_chems <- left_join(ear_sum, select(top_all, site, date, endPoint_top, contr_chems, top), by = c('site', 'date', 'endPoint' = 'endPoint_top')) %>%
    filter(!is.na(top))
  
  plot_summary <- top_all_chems %>%
    group_by(chnm, CAS, Class, contr_chems) %>%
    summarize(median_ear = median(EAR)) %>%
    left_join(top, by = 'contr_chems') %>%
    left_join(select(chem_master, CAS, parent_pesticide, Class, MlWt))
  
  plot_summary$contr_chems <- factor(plot_summary$contr_chems, levels = top$contr_chems)
  
  totals <- plot_summary %>%
    group_by(contr_chems, mix_n_hits_samples) %>%
    summarize(total = sum(median_ear)) %>%
    mutate(chnm = NA)
  
  plot_summary$chnm <- gsub('2,4-Dichlorophenoxyacetic acid', '2,4-D', plot_summary$chnm)
  plot_summary <- ungroup(plot_summary) %>%
    mutate(chnm = ifelse(CAS == '19988-24-0', 'Deethylhydroxyatrazine', chnm)) %>%
    mutate(chnm = ifelse(CAS == '2163-68-0', 'Hydroxyatrazine', chnm)) %>%
    mutate(chnm = ifelse(CAS == '3567-62-2', 'Monomethyldiuron', chnm)) %>%
    mutate(chnm = ifelse(grepl('Deg', Class), paste0('*', chnm), chnm))
  
  # chem order, by parent then by molecular weight
  chem_order <- select(ungroup(plot_summary), parent_pesticide, MlWt, Class, chnm) %>%
    mutate(Class = gsub('Deg - ', '', Class)) %>%
    distinct() %>%
    arrange(Class, parent_pesticide, -MlWt) %>% pull(chnm)
  
  # chem colors
  # adjust alpha for degradates
  colors <- rev(c(RColorBrewer::brewer.pal(8, 'Dark2'), RColorBrewer::brewer.pal(3, 'Set1')[1]))

  colors2 <- c(colors[1:5], colors[6], ggplot2::alpha(colors[6], c(.8, .6, .4, .2)), colors[7], colors[8], ggplot2::alpha(colors[8], c(0.7, 0.4, 0.1)), colors[9])
  plot_summary$chnm <- factor(plot_summary$chnm, levels = chem_order)
  plot_summary$is_parent <- ifelse(plot_summary$chnm == plot_summary$parent_pesticide, 'yes', 'no')

  colors3 <- colors2
  colors3[!grepl('\\*', levels(plot_summary$chnm))] <- 'black'
  colors3[grepl('\\*', levels(plot_summary$chnm))] <- NA
  
  p2 <- ggplot(plot_summary, aes(x = contr_chems, y = median_ear)) +
    geom_bar(stat = 'identity', aes(fill = factor(chnm), color = factor(chnm))) +
    geom_text(data = totals, aes(x = contr_chems, y = total, label = mix_n_hits_samples), 
              vjust = -0.25, color = 'red3', size = 3) +
    theme_bw() +
    geom_hline(yintercept = 0.001, linetype = 2) +
    scale_fill_manual(values = colors2) +
    scale_color_manual(values = colors3) +
    #scale_size_manual(values = c(.22,.6), guide = FALSE) +
    coord_cartesian(ylim = c(0.0002, 0.008)) +
    theme(axis.text.x = element_blank(),
          legend.position = c(0.5, 0.83),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 10)) +
    labs(x = 'Unique mixtures', y = 'Median EARmix', fill = '', color = '') +
    annotate('text', x = 9.5, y = 0.0042, label = 'Number of \nsample hits', color = 'red3', size = 3)

  ggsave(out_file, cowplot::plot_grid(p1, p2), height = 3.6, width = 8)
    
  }

# sites - calculate max EARmix across samples, as well as the 
# number of months where there is an EARmix > 0.001 
calc_site_mix_metrics <- function(top_mixtures, all_samples, EARsum_endpoint, sites_info) {
  
  my_sites <- select(sites_info, site = site_no, shortName)
  
  all_dates <- all_samples %>%
    select(site = SiteID, date = `Sample Date`) %>%
    distinct() %>%
    group_by(site) %>%
    summarize(n_samples = n())
  
  all_possible <- all_samples %>%
    select(site = SiteID, date = `Sample Date`) %>%
    mutate(date = as.POSIXct(date)) %>%
    distinct()
   
  site_mix1 <-  top_mixtures %>%
    group_by(site, date) %>%
    summarize(hits_per_sample = n()) %>%
    right_join(all_possible) %>%
    mutate(hits_per_sample = ifelse(is.na(hits_per_sample), 0, hits_per_sample)) %>%
    group_by(site) %>%
    summarize(med_mix_hits_per_sample = round(median(hits_per_sample), 1),
              mix_hit_n_months = length(unique(lubridate::month(date[hits_per_sample >0]))))
  
  site_mix2 <- EARsum_endpoint %>%
    group_by(site, date) %>%
    summarize(max_sum_ear_endpoint = max(sum_ear_endpoint)) %>%
    right_join(all_possible) %>%
    mutate(max_sum_ear_endpoint = ifelse(is.na(max_sum_ear_endpoint), 0, max_sum_ear_endpoint)) %>%
    group_by(site) %>%
    summarize(med_max_EARmix = round(median(max_sum_ear_endpoint), 4),
              max_EARmix = round(max(max_sum_ear_endpoint), 4))
  
  site_mix <- left_join(site_mix1, site_mix2) %>%
    left_join(my_sites)
  
    
  return(site_mix)
}

calc_top_endpoints <- function(top_mixes) {
  top_ends <- top_mixes %>%
    group_by(endPoint) %>%
    summarize(n_sites = length(unique(site)),
              n_hits = n(), 
              max_sumEAR = max(sum_ear_endpoint),
              parents = paste(unique(unlist(strsplit(contr_parents, split = ', '))), collapse = ', '),
              max_contr_chems = max(n_contr_chems)) %>%
    arrange(-n_hits)
  
  # explore some top mixes
  sub <- filter(top_mixes, endPoint %in% 'NVS_ENZ_hPDE4A1')
  unique(sub$contr_parents)
}
