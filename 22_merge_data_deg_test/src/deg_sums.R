get_parent_class <- function(deg_parent_conc, chem_info) {
  # get class for each parent compound
  parent_class <- select(deg_parent_conc, parent_pesticide, Class, CAS) %>%
    filter(!grepl('Deg', Class)) %>%
    distinct()
  
  other_parents <- chem_info %>%
    filter(!grepl('Deg', Class)) %>%
    filter(!CAS %in% parent_class$CAS) %>%
    select(parent_pesticide = `Chemical Name`, Class, CAS)
  
  all_parents <- bind_rows(parent_class, other_parents)
  
  # some degradates didn't have measured parents
  # so need to manually transfer that information
  
  all_parents <- all_parents %>%
    add_row(parent_pesticide = c('Benomyl', 'Chlorothalonil', 'Parathion-methyl', 'Parathion-ethyl', 'Dacthal'), 
            Class = c(rep('Fungicide',2), rep("Insecticide", 2), 'Herbicide'), 
            CAS = c('17804-35-2', '1897-45-6', '298-00-0', '56-38-2', '1861-32-1'))
  
  return(all_parents)

  }

sum_by_parents <- function(deg_parent_ear, deg_parent_conc, deg_parent_bench) {
  
  #deg_parent_ear$parent_pesticide[is.na(deg_parent_ear$parent_pesticide)] <- as.character(deg_parent_ear$chnm[is.na(deg_parent_ear$parent_pesticide)])
  
 
    
  parent_cmpd_count <- deg_parent_ear %>%
    select(parent_pesticide, chnm) %>%
    distinct() %>%
    group_by(parent_pesticide) %>%
    summarize(n_cmpds = n())
  
 
  
  sum_deg_parents <- function(in_dat, measure) {
    
    # need to account for the different way "benchmark" values are selected
    # if multiple benchmarks per chemical, we're using the most conservative 
    # measure (e.g., highest toxicity, or lowest benchmark value, or max reported TQ)
    if (measure == 'bench') {
      in_dat <- in_dat %>%
        group_by(site, date, CAS, Class, parent_pesticide) %>%
        summarize(EAR = max(EAR))
    }
    parent_sums <- in_dat %>%
      filter(!grepl('Deg', Class)) %>%
      group_by(site, date, parent_pesticide) %>%
      summarize(sumval = sum(EAR)) %>%
      #mutate(type = 'parent') %>%
      ungroup()
    
    deg_sums <- in_dat %>%
      filter(grepl('Deg', Class)) %>%
      group_by(site, date, parent_pesticide) %>%
      summarize(sumval = sum(EAR)) %>%
      #mutate(type = 'degradates') %>%
      ungroup()
    
    j_par_deg_sums <- full_join(rename(parent_sums, p_sumval = sumval), 
                                rename(deg_sums, d_sumval = sumval)) %>%
      rowwise() %>%
      mutate(p_d_sumval = sum(c(p_sumval, d_sumval), na.rm = TRUE))
    
    sum_by_chem <- j_par_deg_sums %>%
      gather(key = 'type', value = 'sumval', p_sumval, d_sumval, p_d_sumval) %>%
      mutate(measure_type = measure)
    
    return(sum_by_chem)
    
  }
  
  ear_graph_deg_sums <- sum_deg_parents(deg_parent_ear, 'ear')
  conc_graph_deg_sums <- sum_deg_parents(deg_parent_conc, 'conc')
  bench_graph_deg_sums <- sum_deg_parents(deg_parent_bench, 'bench')
  
  graph_deg_sums <- bind_rows(ear_graph_deg_sums, conc_graph_deg_sums, bench_graph_deg_sums)
  
  return(graph_deg_sums)
}

crosswalk_names <- function(conc_dat) {
  names_to_cas <- select(conc_dat, chnm, CAS, parent_pesticide, Class) %>%
    distinct()
}
summarize_parents <- function(graph_deg_sums, classes, classes_fixed, zeros) {
  
  if (!zeros) {
    
    graph_deg_sums <- filter(graph_deg_sums, sumval > 0)
  }
  parent_classes <- classes %>%
    filter(parent_pesticide %in% unique(graph_deg_sums$parent_pesticide))
 
  
  sum_by_chem <- graph_deg_sums %>%
    group_by(parent_pesticide, measure_type, type) %>%
    summarize(median_sumval = median(sumval, na.rm = TRUE),
              quant90 = quantile(sumval, probs = 0.90, na.rm = TRUE),
              quant10 = quantile(sumval, probs = 0.10, na.rm = TRUE),
              n_vals = n()) %>%
    #left_join(select(cross, parent_pesticide, CAS)) %>%
    left_join(distinct(select(classes, parent_pesticide, Class))) %>%
    distinct() %>%
    ungroup()
  
  #hits_by_chem <- filter(graph_deg_sums, type == 'p_d_sumval') %>%
  #  summarize
  
  chem_order <- sum_by_chem %>%
    mutate(type_measure = paste0(type, measure_type)) %>%
    select(-quant90, -quant10, -type, -measure_type) %>%
    tidyr::spread(key = type_measure, value = median_sumval) %>%
    rowwise() %>%
    mutate(max = max(c(p_d_sumvalear, p_sumvalear, d_sumvalear), na.rm = TRUE)) %>%
    mutate(max = ifelse(is.infinite(max), NA, max)) %>%
    group_by(Class) %>%
    mutate(max_class = max(max)) %>%
    arrange(max_class, max) %>%
    pull(parent_pesticide) %>% unique()
  
  sum_by_chem$parent_pesticide <- factor(sum_by_chem$parent_pesticide, levels = chem_order)
  
  sum_by_chem$type <- factor(sum_by_chem$type, levels = c('d_sumval', 'p_sumval', 'p_d_sumval'))
  levels(sum_by_chem$type)<- c('degradates', 'parent','parent + degradates')
  
  return(sum_by_chem) 
}

plot_deg_sums <- function(all_dat, top_parents, out_file) {
  
  #all_dat$parent_pesticide[is.na(all_dat$Class)]
  top_dat <- filter(all_dat, parent_pesticide %in% top_parents$parent_pesticide)
  
  line_dat <- data.frame(yintercept = c(0.1, NA, 0.001),
                         measure_type = c('TQchem', 'Conc', 'EARchem'))
  
  top_dat$measure_type <- factor(top_dat$measure_type, levels = c('ear', 'bench', 'conc'))
  levels(top_dat$measure_type) <- c('EARchem', 'TQchem', 'Conc')
  top_dat$Class <- factor(top_dat$Class, levels = c('Herbicide', "Insecticide", "Fungicide", "Other"))
  p <- ggplot(top_dat, aes(y = median_sumval, x = parent_pesticide)) +
    geom_point(aes(color = type), size = 2, alpha = 0.5) +
    facet_grid(rows = vars(Class), cols = vars(measure_type), 
               scales = 'free', shrink = FALSE, space = 'free_y') +
    geom_vline(xintercept = seq(2, length(levels(top_dat$parent_pesticide)), 2), 
               size = 4, color = 'gray88')  +
    geom_point(aes(color = type), size = 2, alpha = 0.5) +
    geom_errorbar(aes(ymin = quant10, ymax = quant90, color = type)) +
    scale_y_log10() +
    coord_flip() +
    geom_hline(data = line_dat, aes(yintercept = yintercept), color = 'orange') +
    theme_bw() +
    labs(x = '', y = '', color = '') +
    theme(legend.position = 'top', legend.direction = 'horizontal')
  
  ggsave(out_file, p, height = 6, width = 6)
  
  
  # ggplot(top_dat, aes(y = median_sumval, x = parent_pesticide)) +
  #   #geom_point(aes(group = type), size = 2, shape = '|',
  #   #           alpha = 0.8, position = position_dodge(width = 1.2)) +
  #   #geom_jitter(aes(color = type), size = 2, alpha = 0.5) +
  #   geom_crossbar(aes(ymin = quant5, ymax = quant95, fill = type), 
  #                 position = position_dodge(width = 1.2), 
  #                 alpha = 0.7) +
  #   facet_wrap(~measure_type) +
  #   #scale_fill_manual(values = c('yellow', 'green', 'blue')) +
  #   scale_y_log10() +
  #   coord_flip() +
  #   geom_hline(yintercept = 0.001)
  # 
  # ggplot(graph_deg_sums, aes(y = median_sumval, x = parent_pesticide)) +
  #   geom_point(aes(group = type, color = type), size = 2, shape = 16,
  #              alpha = 0.8, position = position_dodge(width = 1.2)) +
  #   facet_wrap(~measure_type) +
  #   coord_flip() +
  #   geom_vline(xintercept = seq(2, length(levels(graph_deg_sums$parent_pesticide)), 2), 
  #              size = 4, color = 'gray88') +
  #   geom_point(aes(group = type, color = type), size = 2, shape = 16,
  #              alpha = 0.8, position = position_dodge(width = 1.2)) +
  #   geom_errorbar(aes(ymin = quant5, ymax = quant95, color = type), 
  #                 alpha = 0.8, position = position_dodge(width = 1.2), size = 1) +
  #   #geom_jitter(aes(color = type), size = 2, alpha = 0.5) +
  #   #geom_crossbar(aes(ymin = quant5, ymax = quant95, fill = type), 
  #   #              position = position_dodge(width = 1.2), 
  #   #              alpha = 0.7) +
  #   #scale_fill_manual(values = c('yellow', 'green', 'blue')) +
  #   scale_y_log10() +
  #   
  #   geom_hline(yintercept = 0.001)
  # 
  # ggplot(j_par_deg_sums, aes(y = sumEAR, x = parent_pesticide)) +
  #   geom_boxplot(aes(fill = type)) +
  #   scale_y_log10() +
  #   coord_flip()
}
  

