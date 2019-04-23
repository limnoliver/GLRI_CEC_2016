summarize_deg_parents <- function(deg_parent_ear) {
  
  #deg_parent_ear$parent_pesticide[is.na(deg_parent_ear$parent_pesticide)] <- as.character(deg_parent_ear$chnm[is.na(deg_parent_ear$parent_pesticide)])
  
  parent_class <- select(deg_parent_ear, parent_pesticide, Class) %>%
    filter(!grepl('Deg', Class)) %>%
    distinct()
    
  parent_cmpd_count <- deg_parent_ear %>%
    select(parent_pesticide, chnm) %>%
    distinct() %>%
    group_by(parent_pesticide) %>%
    summarize(n_cmpds = n())
  
  parent_sums <- deg_parent_ear %>%
    filter(!grepl('Deg', Class)) %>%
    group_by(site, date, parent_pesticide) %>%
    summarize(sumEAR = sum(EAR)) %>%
    #mutate(type = 'parent') %>%
    ungroup()
  
  deg_sums <- deg_parent_ear %>%
    filter(grepl('Deg', Class)) %>%
    group_by(site, date, parent_pesticide) %>%
    summarize(sumEAR = sum(EAR)) %>%
    #mutate(type = 'degradates') %>%
    ungroup()
  
  j_par_deg_sums <- full_join(rename(parent_sums, p_sumEAR = sumEAR), 
                              rename(deg_sums, d_sumEAR = sumEAR)) %>%
    rowwise() %>%
    mutate(p_d_sumEAR = sum(c(p_sumEAR, d_sumEAR)))
  
  
  graph_deg_sums <- j_par_deg_sums %>%
    gather(key = 'type', value = 'sumEAR', p_sumEAR, d_sumEAR, p_d_sumEAR) %>%
    group_by(parent_pesticide, type) %>%
    summarize(median_sumEAR = median(sumEAR, na.rm = TRUE),
              quant95 = quantile(sumEAR, probs = 0.95, na.rm = TRUE),
              quant5 = quantile(sumEAR, probs = 0.05, na.rm = TRUE)) %>%
    left_join(parent_class)
  
  chem_order <- graph_deg_sums %>%
    select(-quant95, -quant5) %>%
    tidyr::spread(key = type, value = median_sumEAR) %>%
    rowwise() %>%
    mutate(max = max(c(p_d_sumEAR, p_sumEAR, d_sumEAR), na.rm = TRUE)) %>%
    group_by(Class) %>%
    mutate(max_class = max(max)) %>%
    arrange(max_class, max) %>%
    pull(parent_pesticide)
  
  graph_deg_sums$parent_pesticide <- factor(graph_deg_sums$parent_pesticide, levels = chem_order)
  
  graph_deg_sums$type <- factor(graph_deg_sums$type, levels = c('d_sumEAR', 'p_sumEAR', 'p_d_sumEAR'))
  levels(graph_deg_sums$type)<- c('degradates', 'parent','parent + degradates')
  
  parent_deg_sums <- deg_parent_ear %>%
    left_join(parent_cmpd_count) %>%
    filter(n_cmpds > 1) %>%
    group_by(site, date, parent_pesticide) %>%
    summarize(sumEAR = sum(EAR)) %>%
    mutate(type = 'parent + degradates')
  
  all_dat <- bind_rows(parent_sums, deg_sums, parent_deg_sums) %>%
    left_join(parent_class) %>% 
    group_by(parent_pesticide, type) %>%
    summarize(median_sumEAR = median(sumEAR),
              quant95 = quantile(sumEAR, probs = 0.95),
              quant5 = quantile(sumEAR, probs = 0.05))
  
  test <- filter(all_dat, parent_pesticide %in% 'Tebuthiuron')
  
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
