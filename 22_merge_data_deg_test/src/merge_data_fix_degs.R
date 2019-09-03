find_missing_tox <- function(data_file, parents, metolachlor){
  
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(unique(tox_list$chem_info$CAS))
  ACClong <- remove_flags(ACClong)
  
  # find chems that are not in toxCast
  missing <- unique(tox_list$chem_data$CAS)[-which(unique(tox_list$chem_data$CAS) %in% unique(ToxCast_ACC$CAS))]
  
  return(missing)
}

get_chem_sum_deg <- function(data_file, parents, metolachlor){
  
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(unique(tox_list$chem_info$CAS))
  ACClong <- remove_flags(ACClong)
  
  # find chems that are not in toxCast
  missing <- unique(tox_list$chem_data$CAS)[-which(unique(tox_list$chem_data$CAS) %in% unique(ToxCast_ACC$CAS))]
  
  missing <- filter(tox_list$chem_info, CAS %in% missing)
  
  # crosswalk between missing chems and parent compounds
  missing_parents <- left_join(missing, select(parents, CAS, MlWt, parent_pesticide))

  
  # find out how many sites/samples have missing chem
  missing_sample_counts <- tox_list$chem_data %>%
    filter(CAS %in% missing$CAS) %>%
    group_by(CAS) %>%
    summarize(sample_count = n(),
              site_count = length(unique(SiteID)))
  
 
  
  if (metolachlor == TRUE) {
    missing_parents$parent_pesticide[missing_parents$parent_pesticide == 'Acetochlor/Metolachlor '] <- 'Metolachlor'
  } else {
    missing_parents$parent_pesticide[missing_parents$parent_pesticide == 'Acetochlor/Metolachlor '] <- 'Acetochlor'
    
  }
  
  fixed_deg <- data.frame()
  for (i in 1:nrow(missing_parents)) {
    temp_parent <- as.character(missing_parents$parent_pesticide[i])
    
    if (is.na(temp_parent)) {next}
    
    if (!(temp_parent %in% ACClong$chnm)) {next}
    
    replace_data <- filter(ACClong, chnm %in% temp_parent) %>%
      mutate(CAS = missing_parents$CAS[i],
             chnm = missing_parents$`Chemical Name`[i],
             # account for change on MlWt. Acc is calculated as ACC*MolWt of parent,
             # so need to divide by parent and multiply by degradate MlWt
             ACC_value = ACC_value*(missing_parents$MlWt[i]/MlWt)) %>%
      mutate(MlWt = missing_parents$MlWt[i])
             
    
    fixed_deg <- bind_rows(fixed_deg, replace_data)
  }
  
  ACClong <- bind_rows(ACClong, fixed_deg)
  
  
  cleaned_ep <- clean_endPoint_info(toxEval::end_point_info)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
  
  chemicalSummary <- left_join(chemicalSummary, select(parents, CAS, parent_pesticide)) %>%
    mutate(parent_pesticide = as.character(parent_pesticide))
  
  chemicalSummary$parent_pesticide[is.na(chemicalSummary$parent_pesticide)] <- as.character(chemicalSummary$chnm[is.na(chemicalSummary$parent_pesticide)])
  
  return(chemicalSummary)
}
find_missing_bench <- function(data_file, parents, metolachlor){
  
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  bench_vals <- unique(tox_list$benchmarks) %>%
    filter(!is.na(ACC_value)) %>%
    left_join(select(parents, CAS, MlWt))
  
  # find chems that are not in benchmarks
  missing <- unique(tox_list$chem_data$CAS)[-which(unique(tox_list$chem_data$CAS) %in% unique(bench_vals$CAS))]
  return(missing)
}
get_bench_sum_deg <- function(data_file, parents, metolachlor){
  
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  bench_vals <- unique(tox_list$benchmarks) %>%
    filter(!is.na(ACC_value)) %>%
    left_join(select(parents, CAS, MlWt))
  
  # find chems that are not in benchmarks
  missing <- unique(tox_list$chem_data$CAS)[-which(unique(tox_list$chem_data$CAS) %in% unique(bench_vals$CAS))]
  missing <- filter(tox_list$chem_info, CAS %in% missing)
  
  # crosswalk between missing chems and parent compounds
  missing_parents <- left_join(missing, select(parents, CAS, MlWt, parent_pesticide))
  
  
  # find out how many sites/samples have missing chem
  missing_sample_counts <- tox_list$chem_data %>%
    filter(CAS %in% missing$CAS) %>%
    group_by(CAS) %>%
    summarize(sample_count = n(),
              site_count = length(unique(SiteID)))
  
  
  
  if (metolachlor == TRUE) {
    missing_parents$parent_pesticide[missing_parents$parent_pesticide == 'Acetochlor/Metolachlor '] <- 'Metolachlor'
  } else {
    missing_parents$parent_pesticide[missing_parents$parent_pesticide == 'Acetochlor/Metolachlor '] <- 'Acetochlor'
    
  }
  
  fixed_deg <- data.frame()
  for (i in 1:nrow(missing_parents)) {
    temp_parent <- as.character(missing_parents$parent_pesticide[i])
    
    if (is.na(temp_parent)) {next}
    
    if (!(temp_parent %in% bench_vals$chnm)) {next}
    
    replace_data <- filter(bench_vals, chnm %in% temp_parent) %>%
      mutate(CAS = missing_parents$CAS[i],
             chnm = missing_parents$`Chemical Name`[i],
             # account for change on MlWt. Acc is calculated as ACC*MolWt of parent,
             # so need to divide by parent and multiply by degradate MlWt
             ACC_value = ACC_value*(missing_parents$MlWt[i]/MlWt)) %>%
      mutate(MlWt = missing_parents$MlWt[i])
    
    
    fixed_deg <- bind_rows(fixed_deg, replace_data)
  }
  
  bench_vals <- bind_rows(bench_vals, fixed_deg)
  tox_list$benchmarks <- bench_vals
  
  #cleaned_ep <- clean_endPoint_info(toxEval::end_point_info)
  #filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(tox_list)
  
  chemicalSummary <- left_join(chemicalSummary, select(parents, CAS, parent_pesticide)) %>%
    mutate(parent_pesticide = as.character(parent_pesticide))
  
  chemicalSummary$parent_pesticide[is.na(chemicalSummary$parent_pesticide)] <- as.character(chemicalSummary$chnm[is.na(chemicalSummary$parent_pesticide)])
  
  return(chemicalSummary)
}

get_unclassified <- function(data_file, parents){
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(unique(tox_list$chem_info$CAS))
  ACClong <- remove_flags(ACClong)
  
  # find chems that are not in toxCast
  missing <- unique(tox_list$chem_data$CAS)[-which(unique(tox_list$chem_data$CAS) %in% unique(ToxCast_ACC$CAS))]
  missing <- filter(tox_list$chem_info, CAS %in% missing)
  
  # crosswalk between missing chems and parent compounds
  missing_parents <- left_join(missing, select(parents, CAS, MlWt, parent_pesticide)) %>%
    filter(grepl('/', parent_pesticide))
  
  return(missing_parents)
}

merge_deg_parents <- function(conc_dat, parents, metolachlor) {
  chems <- unique(conc_dat$CAS)
  
  chems_parents <- filter(parents, CAS %in% chems)
  
  conc_dat_parent <- left_join(conc_dat, select(chems_parents, CAS, parent_pesticide)) %>%
    mutate(parent_pesticide = as.character(parent_pesticide))
  
  conc_dat_parent$parent_pesticide[conc_dat_parent$parent_pesticide %in% 'Acetochlor/Metolachlor'] <- ifelse(metolachlor, 'Metolachlor', 'Acetochlor')

  conc_dat_parent$parent_pesticide[is.na(conc_dat_parent$parent_pesticide)] <- as.character(conc_dat_parent$chnm[is.na(conc_dat_parent$parent_pesticide)])
  
  
  
  return(conc_dat_parent)
}

plot_parent_deg <- function(conc_dat) {
  meto <- select(conc_dat, -CAS) %>%
    filter(parent_pesticide == 'Metolachlor') %>%
    tidyr::spread(key = chnm, value = EAR)
  
  ggplot(meto, aes(x = ))
  ggplot(meto, aes(x = shortName, y = EAR)) +
    geom_boxplot(aes(fill = chnm)) +
    scale_y_log10()
      
  
  library(ggplot2)
  
  head(conc_dat)
  par_deg <- conc_dat %>%
    select(chnm, parent_pesticide) %>%
    distinct() %>%
    group_by(parent_pesticide) %>%
    summarize(n_compounds = n()) %>%
    filter(n_compounds > 1)
  
  plot_par_deg <- filter(conc_dat, parent_pesticide %in% par_deg$parent_pesticide) %>%
    filter(!is.na(parent_pesticide)) %>%
    mutate(parent = ifelse(grepl('deg', Class, ignore.case = TRUE), FALSE, TRUE))
                
  
  ggplot(plot_par_deg, aes(x = as.factor(chnm), y = EAR)) +
    geom_boxplot(aes(fill = Class), position = position_dodge2(preserve = "total"), 
                 varwidth = FALSE, width = 0.4) +
    facet_wrap(~parent_pesticide, ncol = 1, drop = TRUE, scales = 'free_y', shrink = FALSE, dir = 'v') +
    coord_flip() +
    scale_y_log10() +
    theme(strip.text.x = element_blank(),
          axis.text.y = element_text(size = 6))
}

plot_unclassified_degs <- function(conc_dat, unclassified_degs, out_file) {
  metolachlor_cas <- filter(conc_dat, chnm %in% 'Metolachlor') %>% 
    select(CAS) %>%
    distinct()
  
  acetochlor_cas <- filter(conc_dat, chnm %in% 'Acetochlor') %>% 
    select(CAS) %>%
    distinct()
  
  deg <- filter(conc_dat, CAS %in% unclassified_degs$CAS)
  met <- filter(conc_dat, CAS %in% metolachlor_cas$CAS) %>%
    rename(met_ear = EAR)
  acet <- filter(conc_dat, CAS %in% acetochlor_cas$CAS ) %>%
    rename(acet_ear = EAR)
  
  deg_parent <- left_join(deg, select(met, site, date, met_ear)) %>%
    left_join(select(acet, site, date, acet_ear))
  
  p1 <- ggplot(deg_parent, aes(x = acet_ear, y = EAR)) +
    geom_point(aes(color = chnm)) +
    theme_bw() +
    labs(x = 'Acetochlor concentration', y = 'Degradate concentration')
  
  p2 <- ggplot(deg_parent, aes(x = met_ear, y = EAR)) +
    geom_point(aes(color = chnm)) +
    theme_bw() +
    labs(x = 'Metolochlor concentration', y = 'Degradate concentration')
  
  
  ggsave(out_file, cowplot::plot_grid(p1, p2, align = 'h', nrow = 2), height = 6, width = 7)
    
}