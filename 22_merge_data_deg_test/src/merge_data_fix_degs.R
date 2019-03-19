get_chem_sum_deg <- function(data_file, parents, metolachlor){
  
  tox_list <- create_toxEval(data_file)
  tox_list$chem_data <- filter(tox_list$chem_data, Value != 0)
  
  ACClong <- get_ACC(tox_list$chem_info$CAS)
  ACClong <- remove_flags(ACClong)
  
  # find chems that are not in toxCast
  missing <- unique(tox_list$chem_data$CAS)[-which(unique(tox_list$chem_data$CAS) %in% unique(ToxCast_ACC$CAS))]
  missing <- filter(tox_list$chem_info, CAS %in% missing)
  
  # crosswalk between missing chems and parent compounds
  missing_parents <- left_join(missing, select(parents, CAS, MlWt, parent_pesticide))
  
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
             MlWt = missing_parents$MlWt[i])
    
    fixed_deg <- bind_rows(fixed_deg, replace_data)
  }
  
  ACClong <- bind_rows(ACClong, fixed_deg)
  
  
  cleaned_ep <- clean_endPoint_info(toxEval::end_point_info)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
  return(chemicalSummary)
}

merge_deg_parents <- function(conc_dat, parents, metolachlor) {
  chems <- unique(conc_dat$CAS)
  
  chems_parents <- filter(parents, CAS %in% chems)
  
  chems_parents$parent_pesticide[chems_parents$parent_pesticide == 'Acetochlor/Metolachlor'] <- ifelse(metolachlor, 'Metolachlor', 'Acetochlor')

  conc_dat_parent <- left_join(conc_dat, select(parents, CAS, parent_pesticide))
  
  return(conc_dat_parent)
}

plot_parent_deg <- function(conc_dat) {
  head(conc_dat)
  
  library(ggplot2)
  
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
          axis.text.x = element_text(size = 8))
}