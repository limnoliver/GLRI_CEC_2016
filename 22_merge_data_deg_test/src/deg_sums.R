get_parent_class <- function(deg_parent_conc) {
  # get class for each parent compound
  parent_class <- select(deg_parent_conc, parent_pesticide, Class) %>%
    filter(!grepl('Deg', Class)) %>%
    distinct()
}

sum_by_parents <- function(deg_parent_ear, deg_parent_conc, deg_parent_bench) {
  
  #deg_parent_ear$parent_pesticide[is.na(deg_parent_ear$parent_pesticide)] <- as.character(deg_parent_ear$chnm[is.na(deg_parent_ear$parent_pesticide)])
  
 
    
  parent_cmpd_count <- deg_parent_ear %>%
    select(parent_pesticide, chnm) %>%
    distinct() %>%
    group_by(parent_pesticide) %>%
    summarize(n_cmpds = n())
  
 
  
  sum_deg_parents <- function(in_dat, measure) {
    
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
      mutate(p_d_sumval = sum(c(p_sumval, d_sumval)))
    
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

summarize_parents <- function(graph_deg_sums, classes) {
  
  sum_by_chem <- graph_deg_sums %>%
    group_by(parent_pesticide, measure_type, type) %>%
    summarize(median_sumval = median(sumval, na.rm = TRUE),
              quant95 = quantile(sumval, probs = 0.95, na.rm = TRUE),
              quant5 = quantile(sumval, probs = 0.05, na.rm = TRUE)) %>%
    left_join(classes)
  
  chem_order <- sum_by_chem %>%
    mutate(type_measure = paste0(type, measure_type)) %>%
    select(-quant95, -quant5, -type, -measure_type) %>%
    tidyr::spread(key = type_measure, value = median_sumval) %>%
    rowwise() %>%
    mutate(max = max(c(p_d_sumvalear, p_sumvalear, d_sumvalear), na.rm = TRUE)) %>%
    group_by(Class) %>%
    mutate(max_class = max(max)) %>%
    arrange(max_class, max) %>%
    pull(parent_pesticide)
  
  graph_deg_sums$parent_pesticide <- factor(graph_deg_sums$parent_pesticide, levels = chem_order)
  
  graph_deg_sums$type <- factor(graph_deg_sums$type, levels = c('d_sumEAR', 'p_sumEAR', 'p_d_sumEAR'))
  levels(graph_deg_sums$type)<- c('degradates', 'parent','parent + degradates')
  
  return(graph_deg_sums) 
}

plot_deg_sums <- function() {
  
  library(ggplot2)
  ggplot(all_dat, aes(y = median_sumEAR, x = parent_pesticide)) +
    geom_point(aes(color = type), size = 2, alpha = 0.5) +
    geom_errorbar(aes(ymin = quant5, ymax = quant95, color = type)) +
    scale_y_log10() +
    coord_flip()
  
  ggplot(graph_deg_sums, aes(y = median_sumEAR, x = parent_pesticide)) +
    #geom_point(aes(group = type), size = 2, shape = '|',
    #           alpha = 0.8, position = position_dodge(width = 1.2)) +
    #geom_jitter(aes(color = type), size = 2, alpha = 0.5) +
    geom_crossbar(aes(ymin = quant5, ymax = quant95, fill = type), 
                  position = position_dodge(width = 1.2), 
                  alpha = 0.7) +
    #scale_fill_manual(values = c('yellow', 'green', 'blue')) +
    scale_y_log10() +
    coord_flip() +
    geom_hline(yintercept = 0.001)
  
  ggplot(graph_deg_sums, aes(y = median_sumEAR, x = parent_pesticide)) +
    geom_point(aes(group = type, color = type), size = 2, shape = 16,
               alpha = 0.8, position = position_dodge(width = 1.2)) +
    
    coord_flip() +
    geom_vline(xintercept = seq(2, length(levels(graph_deg_sums$parent_pesticide)), 2), 
               size = 4, color = 'gray88') +
    geom_point(aes(group = type, color = type), size = 2, shape = 16,
               alpha = 0.8, position = position_dodge(width = 1.2)) +
    geom_errorbar(aes(ymin = quant5, ymax = quant95, color = type), 
                  alpha = 0.8, position = position_dodge(width = 1.2), size = 1) +
    #geom_jitter(aes(color = type), size = 2, alpha = 0.5) +
    #geom_crossbar(aes(ymin = quant5, ymax = quant95, fill = type), 
    #              position = position_dodge(width = 1.2), 
    #              alpha = 0.7) +
    #scale_fill_manual(values = c('yellow', 'green', 'blue')) +
    scale_y_log10() +
    
    geom_hline(yintercept = 0.001)
  
  ggplot(j_par_deg_sums, aes(y = sumEAR, x = parent_pesticide)) +
    geom_boxplot(aes(fill = type)) +
    scale_y_log10() +
    coord_flip()
}
  

