find_missing_tox <- function(chem_conc, chem_ear, chem_info){
  
  detected <- unique(chem_conc$CAS)
  in_toxcast <- unique(chem_ear$CAS)
  
  
  chems_missing <- chem_info %>%
    filter(CAS %in% detected) %>%
    filter(!CAS %in% in_toxcast)
  
  return(chems_missing)
}

complete_parents <- function(chem_master) {
  # some parents are not in this study so don't have additional information
  # this function adds complete parent info for all compounds
  parents <- filter(chem_master, `Chemical Name` == parent_pesticide | compound == parent_pesticide)
  
  # now find compounds with no parent rep
  lost_parents <- filter(chem_master, !(parent_pesticide %in% parents$`Chemical Name` | parent_pesticide %in% parents$compound))
  
  parents <- parents %>%
    select(parent_pesticide, parent_CAS = CAS, parent_MlWt = MlWt)
    
  
  return(parents)
}

get_chem_sum_deg <- function(data_file, missing_chems, parents, metolachlor, chem_master){
 
  
  if (metolachlor == TRUE) {
    missing_chems$parent_pesticide[missing_chems$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Metolachlor'
  } else {
    missing_chems$parent_pesticide[missing_chems$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Acetochlor'
    
  }
  
  
  # some missing compounds are parent compounds, and are not applicable to this exercise
  #missing_chems_degs
  
  # three chemicals were not matched because parent was not measured in this study
  # see if the parent is in toxCast
  # all three are. Maybe put this content elsewhere?
  missing_chem_parents <- left_join(missing_chems, parents) %>%
    filter(!(`Chemical Name` == parent_pesticide | compound == parent_pesticide))

  ACClong_parents <- get_ACC(unique(missing_chem_parents$parent_CAS))
  ACClong_parents <- remove_flags(ACClong_parents, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical", "ACCLessThan"))
  
  
  fixed_deg <- data.frame()
  
  for (i in 1:nrow(missing_chem_parents)) {

    if (is.na(missing_chem_parents$parent_CAS[i])) {next}
    
    if (!(missing_chem_parents$parent_CAS[i] %in% ACClong_parents$CAS)) {next}
    
    replace_data <- filter(ACClong_parents, CAS %in% missing_chem_parents$parent_CAS[i]) %>%
      mutate(CAS = missing_chem_parents$CAS[i],
             chnm = missing_chem_parents$`Chemical Name`[i],
             # account for change on MlWt. Acc is calculated as ACC*MolWt of parent,
             # so need to divide by parent and multiply by degradate MlWt
             ACC_value = ACC_value*(missing_chem_parents$MlWt[i]/MlWt)) %>%
      mutate(MlWt = missing_chem_parents$MlWt[i])
             
    
    fixed_deg <- bind_rows(fixed_deg, replace_data)
  }
  
  
  # ACC vals for all other compounds
  tox_list <- create_toxEval(data_file)
  ACClong <- get_ACC(unique(tox_list$chem_info$CAS))
  ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical", 'ACCLessThan'))
  
  # have to carry exclusions for the degradates through, 
  # given parent exclusions
  exclude <- tox_list$exclusions %>%
    filter(CAS %in% missing_chem_parents$parent_CAS) %>%
    rename(parent_CAS = CAS)
  
  missing_exclude <- select(missing_chem_parents, CAS, parent_CAS) %>%
    filter(parent_CAS %in% unique(exclude$parent_CAS)) %>%
    left_join(exclude) %>% select(-parent_CAS)
  
  tox_list$exclusions <- bind_rows(tox_list$exclusions, missing_exclude)
    
    
  
  
  # verify there is no overlap between replacement data and original data
  all(!unique(fixed_deg$CAS) %in% unique(ACClong$CAS))
  
  ACClong <- bind_rows(ACClong, fixed_deg)
  
  
  cleaned_ep <- clean_endPoint_info(toxEval::end_point_info)
  filtered_ep <- filter_groups(cleaned_ep)

  chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
  
  
  chemicalSummary <- left_join(chemicalSummary, select(chem_master, CAS, parent_pesticide))
  
  if (metolachlor == TRUE) {
    chemicalSummary$parent_pesticide[chemicalSummary$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Metolachlor'
  } else {
    chemicalSummary$parent_pesticide[chemicalSummary$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Acetochlor'
    
  }

  return(chemicalSummary)
}
find_missing_bench <- function(chem_conc, chem_bench, chem_info){
  detected <- unique(chem_conc$CAS)
  in_bench <- unique(chem_bench$CAS)
  
  length(in_bench)
  
  chems_missing <- chem_info %>%
    filter(CAS %in% detected) %>%
    filter(!CAS %in% in_bench)

  return(chems_missing)
}
get_bench_sum_deg <- function(data_file, missing_chems, parents, metolachlor, chem_master){
  
  if (metolachlor == TRUE) {
    missing_chems$parent_pesticide[missing_chems$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Metolachlor'
  } else {
    missing_chems$parent_pesticide[missing_chems$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Acetochlor'
    
  }
  # reduce missing chems to just deg compounds
  # some missing are parents
  missing_chem_parents <- left_join(missing_chems, parents) %>%
    filter(!(`Chemical Name` == parent_pesticide | compound == parent_pesticide))
  
  tox_list <- create_toxEval(data_file)
  
  bench_vals <- distinct(tox_list$benchmarks) %>%
    filter(!is.na(ACC_value))
  
  # note I manually checked if there were any missing compounds in benchmarks that
  # did not have CAS numbers that matched compounds in our database without CAS numbers
  # I found no matches so did not try to do a name match.
  
  fixed_deg <- data.frame()
  for (i in 1:nrow(missing_chem_parents)) {

    if (is.na(missing_chem_parents$parent_pesticide[i])) {next}
    
    if (!(missing_chem_parents$parent_CAS[i] %in% unique(bench_vals$CAS))) {next}
    
    replace_data <- filter(bench_vals, CAS %in% missing_chem_parents$parent_CAS[i]) %>%
      mutate(CAS = missing_chem_parents$CAS[i],
             chnm = missing_chem_parents$`Chemical Name`[i],
             orig_name = missing_chem_parents$compound[i],
             # account for change on MlWt. Acc is calculated as ACC*MolWt of parent,
             # so need to divide by parent and multiply by degradate MlWt
             ACC_value = ACC_value*(missing_chem_parents$MlWt[i]/missing_chem_parents$parent_MlWt[i])) %>%
      mutate(MlWt = missing_chem_parents$MlWt[i])
    
    
    fixed_deg <- bind_rows(fixed_deg, replace_data)
  }
  
  bench_vals <- bind_rows(bench_vals, fixed_deg)
  tox_list$benchmarks <- bench_vals
  
  #cleaned_ep <- clean_endPoint_info(toxEval::end_point_info)
  #filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(tox_list)
  
  chemicalSummary <- left_join(chemicalSummary, select(chem_master, CAS, parent_pesticide))
  
  if (metolachlor == TRUE) {
    chemicalSummary$parent_pesticide[chemicalSummary$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Metolachlor'
  } else {
    chemicalSummary$parent_pesticide[chemicalSummary$parent_pesticide == 'Acetochlor/Metolachlor'] <- 'Acetochlor'
    
  }
  
  return(chemicalSummary)
}

create_tox_supp_table <- function(toxcast_dat, bench_dat, all_chems) {
  # which compounds are measured in toxcast
  # which compounds had assumed parent vals in toxcast
  
  eligible_chems <- all_chems$CAS[!is.na(all_chems$CAS)]
  
  in_toxcast <- all_chems %>%
    filter(!is.na(CAS)) %>%
    filter(CAS %in% unique(toxEval::ToxCast_ACC$CAS))
  
  
  used_toxcast <- unique(toxcast_dat$CAS)[unique(toxcast_dat$CAS) %in% unique(toxEval::ToxCast_ACC$CAS)]
  used_parent_acc <- unique(toxcast_dat$CAS)[!unique(toxcast_dat$CAS) %in% unique(toxEval::ToxCast_ACC$CAS)]
  
  # which compounds are measured in benchmarks
  # which compounds had assumed parent valus in benchmarks
}

get_unclassified <- function(data_file, parents){
  tox_list <- create_toxEval(data_file)
  #tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(unique(tox_list$chem_info$CAS))
  ACClong <- remove_flags(ACClong, flagsShort = c("Borderline", "OnlyHighest", "GainAC50", "Biochemical", 'ACCLessThan'))
  
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