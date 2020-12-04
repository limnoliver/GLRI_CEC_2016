# get chemicals that are in >50% of months at >50% of sites
library(lubridate)
select_top_chems <- function(chem_dat) {
  count_by_site_month <- mutate(chem_dat, month = month(date)) %>%
    group_by(site, CAS, chnm) %>%
    summarize(month_count = length(unique(month)))
  
  chems_count_sites <- group_by(count_by_site_month, CAS, chnm) %>%
    summarize(n_sites_2months = length(which(month_count >= 2)),
              n_sites = n())
  
  top_chems <- filter(chems_count_sites, n_sites_2months >= 8)
  return(top_chems)
  
  top_chems_attr <- left_join(top_chems, distinct(select(chem_dat, CAS, Class)))
}
calc_chem_tox_hits <- function(chem_tox_dat, chem_conc_dat) {
  
  site_hits_tox <- chem_tox_dat %>%
    group_by(site, date, CAS, chnm) %>%
    summarize(sumEAR = sum(EAR)) %>%
    filter(sumEAR > 0.001) %>%
    ungroup() %>%
    group_by(CAS, chnm) %>%
    summarize(n_hits = n(),
              n_site_hits = length(unique(site))) %>%
    filter(n_site_hits >= sites_hit)
  
  # site_samples <- reduced_dat %>%
  #   group_by(SiteID, sample_dt) %>%
  #   summarize(n_chems = n()) %>%
  #   group_by(SiteID) %>%
  #   summarize(n_samples = n()) %>%
  #   rename(site = SiteID)
  # 
  # hits_per_site <- chem_tox_dat %>%
  #   group_by(site, date, CAS, chnm) %>%
  #   summarize(sumEAR = sum(EAR)) %>%
  #   filter(sumEAR > 0.001) %>%
  #   ungroup() %>%
  #   group_by(site, CAS, chnm) %>%
  #   summarize(n_hits_per_site = n()) %>%
  #   left_join(site_samples) %>%
  #   mutate(prop_sample_hits = n_hits_per_site/n_samples) %>%
  #   group_by(CAS, chnm) %>%
  #   summarize(mean_prop_sample_hits = mean(prop_sample_hits))
    
  
  site_detects <- chem_conc_dat %>%
    group_by(CAS, chnm) %>%
    summarize(n_sites = length(unique(site))) %>%
    filter(n_sites >= sites_detect)
  
  top_chems <- inner_join(site_hits_tox, site_detects, by = c('CAS', 'chnm'))
  return(top_chems)
}

calc_chem_tox_metrics <- function(chem_dat) {
 
  hits_tox <- chem_dat %>%
    group_by(site, date, CAS, chnm) %>%
    summarize(sumEAR = sum(EAR)) %>%
    group_by(CAS, chnm) %>%
    summarize(tox_n_detect = n(), 
              n_hits_tox = length(which(sumEAR > 0.01)), 
              n_sites_tox_hits = length(unique(site[which(sumEAR > 0.01)])),
              mean_sumEAR = mean(sumEAR),
              max_sumEAR = max(sumEAR))
  
  return(hits_tox)
}

calc_chem_wq_metrics <- function(chem_dat) {

  hits_wq <- chem_dat %>%
    group_by(site, date, CAS, chnm) %>%
    summarize(max_bench = max(EAR)) %>%
    group_by(CAS, chnm) %>%
    summarize(wq_n_detect = n(), 
              n_hits_wq = length(which(max_bench > 0.1)), 
              n_sites_hits = length(unique(site[which(max_bench > 0.1)])),
              mean_bench = mean(max_bench),
              max_bench = max(max_bench))
  
  return(hits_wq)
  
}

cal_chem_mix_metrics <- function(mix_dat) {
  
}

merge_top_chems_hits <- function(tox_hits, wq_hits, top_chems) {
  dat <- left_join(top_chems, wq_hits) %>%
    left_join(tox_hits)
  
  return(dat)
}

summarize_chems <- function(file_name, chem_vals, chem_crosswalk, chem_info, conc_dat) {
  sample_count <- all_chems %>%
    group_by(CAS, pCode) %>%
    summarize(n_samples = n(),
              n_sites = length(unique(SiteID)))
  
  n_detect <- all_chems %>%
    filter(!remark_cd %in% '<') %>%
    group_by(CAS) %>%
    summarize(n_detect = n(),
              n_detect_sites = length(unique(SiteID)))
  
  detects <- left_join(sample_count, n_detect) %>%
    mutate(n_detect = ifelse(is.na(n_detect), 0, n_detect),
           n_detect_sites = ifelse(is.na(n_detect_sites), 0, n_detect_sites),
           perc_detect = round(100*(n_detect/n_samples), 0),
           perc_detect_site = round(100*(n_detect_sites/n_sites), 0)) %>%
    mutate(detection_column = paste0(perc_detect, ' (', n_detect, '/', n_samples, ')'),
           detection_site_column = paste0(perc_detect_site, ' (', n_detect_sites, '/', n_sites, ')')) %>%
    select(CAS, pCode, detection_column, detection_site_column) %>%
    left_join(select(chem_crosswalk, CAS, compound, parent_pesticide)) %>%
    left_join(chem_info)
  
  write.csv(detects, file_name, row.names = FALSE)
    
}

summarize_chem_meta <- function(file_name, chem_vals, chem_crosswalk, chem_info_all, chems_missing_toxcast, chems_missing_bench, neonic) {
  
  sample_count <- chem_vals %>%
    group_by(CAS, pCode) %>%
    summarize(n_samples = n(),
              n_sites = length(unique(SiteID)))
  
  n_detect <- chem_vals %>%
    filter(!remark_cd %in% '<') %>%
    group_by(CAS) %>%
    summarize(n_detect = n(),
              n_detect_sites = length(unique(SiteID))) %>%
    mutate(n_detect = ifelse(is.na(n_detect), 0, n_detect),
           n_detect_sites = ifelse(is.na(n_detect_sites), 0, n_detect_sites))

  meta <- left_join(sample_count, chem_info_all, by = 'CAS') %>%
    left_join(select(chem_crosswalk, CAS, compound, parent_pesticide)) %>%
    mutate(`In toxCast` = ifelse(CAS %in% chems_missing_toxcast, '', 'x'),
           `In benchmarks` = ifelse(CAS %in% chems_missing_bench, '', 'x')) %>%
    left_join(n_detect, by = 'CAS')
  
  meta$parent_pesticide[is.na(meta$parent_pesticide)] <- meta$`Chemical Name`[is.na(meta$parent_pesticide)]
  
  meta <- meta %>%
    arrange(parent_pesticide) %>%
    select(`Chemical Name`, CAS, `USGS parameter code` = pCode, 
           Class, `Parent compound` = parent_pesticide, `In toxCast`, 
           `In benchmarks`, `N samples measured` = n_samples, `N samples detected` = n_detect, 
           `N sites measured` = n_sites, `N sites detected` = n_detect_sites)
  
  
  write.csv(meta, file_name, row.names = FALSE)
  
}

calc_parent_tox_hits <- function(parent_sums, mixtures, all_dat, chem_master){
  #all_dat <- make('chem_data_complete')
  # every sample
  all_samples <- all_dat %>%
    left_join(select(chem_master, pCode, parent_pesticide)) %>%
    group_by(parent_pesticide, SiteID, `Sample Date`) %>%
    summarize(n_chems = n()) %>%
    group_by(parent_pesticide) %>%
    summarize(n_samples = n(),
              n_sites  = length(unique(SiteID)))

  sample_dates_ear <- all_dat %>%
    left_join(select(chem_master, pCode, parent_pesticide)) %>%
    mutate(`Sample Date`=as.POSIXct(`Sample Date`)) %>%
    select(site = SiteID, date = `Sample Date`, parent_pesticide) %>%
    filter(parent_pesticide %in% unique(parent_sums$parent_pesticide[parent_sums$measure_type == 'ear'])) %>%
    distinct()
  
  sample_dates_bench <- all_dat %>%
    left_join(select(chem_master, pCode, parent_pesticide)) %>%
    mutate(`Sample Date`=as.POSIXct(`Sample Date`)) %>%
    select(site = SiteID, date = `Sample Date`, parent_pesticide) %>%
    filter(parent_pesticide %in% unique(parent_sums$parent_pesticide[parent_sums$measure_type == 'bench'])) %>%
    distinct()
  
  sample_dates_conc <- all_dat %>%
    left_join(select(chem_master, pCode, parent_pesticide)) %>%
    mutate(`Sample Date`=as.POSIXct(`Sample Date`)) %>%
    select(site = SiteID, date = `Sample Date`, parent_pesticide) %>%
    filter(parent_pesticide %in% unique(parent_sums$parent_pesticide[parent_sums$measure_type == 'conc'])) %>%
    distinct()
  
  
  
  hits_ear <- parent_sums %>%
    filter(type == 'p_d_sumval' & measure_type %in% 'ear') %>%
    right_join(sample_dates_ear) %>%
    left_join(ungroup(all_samples)) %>%
    mutate(type = 'p_d_sumval', measure_type = 'ear',
           sumval = ifelse(is.na(sumval), 0, sumval)) %>%
    #filter(sumval > 0.001) %>%
    group_by(parent_pesticide) %>%
    summarize(ear_hit_prob = round(length(which(sumval > 0.001))/unique(n_samples),2), 
              ear_hit_sites_prob = round(length(unique(site[sumval > 0.001]))/unique(n_sites), 2),
              ear_hit_months = length(unique(lubridate::month(date[sumval > 0.001]))),
              ear_median_sumval = round(median(sumval), 4),
              ear_max_sumval = round(max(sumval), 4)) 
  
  hits_bench <- parent_sums %>%
    filter(type == 'p_d_sumval' & measure_type %in% 'bench') %>%
    right_join(sample_dates_bench) %>%
    left_join(ungroup(all_samples)) %>%
    mutate(type = 'p_d_sumval', measure_type = 'bench',
           sumval = ifelse(is.na(sumval), 0, sumval)) %>%
    #filter(sumval > 0.01) %>%
    group_by(parent_pesticide) %>%
    summarize(bench_hit_prob = round(length(which(sumval > 0.1))/unique(n_samples), 2), 
              bench_hit_sites_prob = round(length(unique(site[sumval > 0.1]))/unique(n_sites), 2),
              bench_hit_months = length(unique(lubridate::month(date[sumval > 0.1]))),
              bench_median_sumval = round(median(sumval), 3),
              bench_max_sumval = round(max(sumval), 3))
  
  hits_mixes <- mixtures
  
  detections <- parent_sums %>%
    left_join(ungroup(all_samples)) %>%
    filter(type == 'p_d_sumval' & measure_type %in% 'conc') %>%
    group_by(parent_pesticide) %>%
    summarize(detect_prob = round(length(which(sumval > 0))/unique(n_samples), 2),
              detect_sites_prob = round(length(unique(site[sumval > 0]))/(unique(n_sites)),2),
              detect_months = length(unique(lubridate::month(date[sumval > 0]))))
  
  parent_metrics <- left_join(detections, hits_mixes, by = 'parent_pesticide') %>%
    left_join(hits_bench, by = 'parent_pesticide') %>%
    left_join(hits_ear, by = 'parent_pesticide')
  
  return(parent_metrics)
    
}

determine_priorities <- function(metrics, 
                                 site_hits, sample_hits, occurrence, 
                                 missing_toxcast, missing_bench, mixtures) {
  
  # top_compounds <- ranks %>%
  #   mutate(priority = case_when(
  #     mean_rel_val >= bench & metric_type %in% 'TQchem' ~ TRUE,
  #     mean_rel_val >= ear & metric_type %in% c('EARmix', "EARchem") ~ TRUE,
  #     mean_rel_val >= occurrence & metric_type %in% 'Occurrence' ~ TRUE,
  #     TRUE ~ FALSE
  #   ))
  metric_long <- metrics %>%
    tidyr::gather(key = 'metric', value = 'value', -parent_pesticide)
  
  return_top <- function(focal_metric, top_value) {
    top_mix <- filter(metric_long, metric %in% focal_metric) %>% 
      filter(value >= top_value) %>% 
      pull(parent_pesticide)
    
    return(top_mix)
  }
  
  top_ear <- unique(c(return_top('ear_hit_prob', sample_hits), return_top('ear_hit_sites_prob', site_hits)))
  top_bench <- unique(c(return_top('bench_hit_prob', sample_hits), return_top('bench_hit_sites_prob', site_hits)))
  
  # now find top based on occurance
  # filter to those missing tox AND bench
  missing <- bind_rows(missing_bench, missing_toxcast) %>%
    group_by(CAS, pCode, `Chemical Name`, parent_pesticide) %>%
    summarize(n_missing = n()) %>% 
    filter(n_missing > 1)

  # which of these missing compounds is also detected at high frequency? 
  browser()
  occur_stats <- occurrence %>%
    filter(CAS %in% unique(missing$CAS)) %>%
    group_by(CAS) %>%
    summarize(n_measured = n(),
              n_detected = sum(!remark_cd %in% '<'),
              n_sites_measured = length(unique(SiteID)),
              n_sites_detected = length(unique(SiteID[!remark_cd %in% '<']))) %>%
    mutate(sites_prob = round(n_sites_detected/n_sites_measured, 2),
           samples_prob = round(n_detected/n_measured, 2)) %>%
    filter(sites_prob >= 2*site_hits | samples_prob >= 2*sample_hits)

 top_occur <- filter(missing, CAS %in% occur_stats$CAS) %>%
    pull(parent_pesticide)
  
  top_mix <- mixtures %>%
    filter(mixes_essential_hit_prob >=sample_hits | mixes_essential_hit_sites_prob >= site_hits) %>%
    pull(parent_pesticide)
  
  top_compounds <- data.frame(parent_pesticide = c(top_ear, top_bench, top_occur, top_mix),
                              metric_type = c(rep('EARchem', length(top_ear)), 
                                              rep('TQchem', length(top_bench)),
                                              rep('Occurrence', length(top_occur)),
                                              rep('EARmix', length(top_mix))))
 

  
  return(top_compounds)
  
}

filter_top <- function(top) {
  #top_names <- filter(top, priority)
  
  return(unique(top$parent_pesticide))
}

rank_parents <- function(parent_metrics) {

  metrics_long <- parent_metrics %>%
    tidyr::gather(key = 'metric', value = 'metric_value', -parent_pesticide) %>%
    mutate(metric_type = case_when(
      grepl('mixes', metric) ~ 'EARmix',
      grepl('detect', metric) ~ 'Occurrence',
      grepl('ear', metric) ~ 'EARchem',
      grepl('bench', metric) ~ 'TQchem')) %>%
    group_by(metric, metric_type) %>%
    mutate(rank_value = rank(-metric_value),
           relative_value = metric_value/max(metric_value, na.rm = TRUE))
  
 return(metrics_long)
  

}

calc_mean_ranks <- function(ranks) {
  metrics_groups <- ungroup(ranks) %>%
    group_by(parent_pesticide, metric_type) %>%
    summarize(mean_rel_val = round(mean(relative_value), 2))
  
  return(metrics_groups)
}

plot_top_parents <- function(top, metadata, conc, ear, bench, raw_ear, raw_bench, out_file) {
  meta <- select(metadata, parent_pesticide, Class) %>%
    distinct() %>%
    filter(!grepl('Deg - ', Class))
  plot_dat <- left_join(top, meta)
  plot_dat$Class <- factor(plot_dat$Class, c('Herbicide', 'Fungicide', 'Insecticide'))
  plot_dat$metric_type <- factor(plot_dat$metric_type, c('EARmix', 'EARchem', 'TQchem', 'Occurrence'))
  
  
  keep_chems <- select(top, parent_pesticide) %>% distinct()
  
  plot_dat <- filter(plot_dat, parent_pesticide %in% keep_chems$parent_pesticide)
  
  # second plot data
  all_samples <- conc %>%
    select(parent_pesticide, site, date) %>%
    distinct()
  
  all_compounds <- conc %>% 
    select(CAS, chnm, parent_pesticide) %>%
    distinct()
  
  all_possible <- full_join(all_samples, all_compounds)
  
  conc_median <- left_join(all_possible, conc) %>%
    mutate(EAR = ifelse(is.na(EAR), 0, EAR)) %>%
    group_by(CAS, chnm, parent_pesticide) %>%
    summarize(median_val = median(EAR), n = n()) %>%
    mutate(metric = 'Concentration',
           estimated = FALSE)

  # was the EAR value estimated? 
  ear_true <- raw_ear %>%
    select(chnm) %>%
    distinct()
  
  ear_sum<- ear %>%
    group_by(chnm, CAS, parent_pesticide, site, date) %>%
    summarize(sumval = sum(EAR))
  
  
  # fill in zeros when a parent/or deg wasn't measured if another was
  all_samples <- ear %>%
    select(parent_pesticide, site, date) %>%
    distinct()
  
  all_compounds <- ear %>% 
    select(CAS, chnm, parent_pesticide) %>%
    distinct()
  
  all_possible <- full_join(all_samples, all_compounds)
  
  ear_median <- left_join(all_possible, ungroup(ear_sum)) %>%
    mutate(sumval = ifelse(is.na(sumval), 0, sumval)) %>%
    group_by(CAS, chnm, parent_pesticide) %>%
    summarize(median_val = median(sumval)) %>%
    mutate(metric = 'EARchem',
           estimated = ifelse(chnm %in% ear_true$chnm, FALSE, TRUE))
  
  # use EAR median to order chems
  
  order_chems <- group_by(ear_median, parent_pesticide) %>%
    summarize(sum_median_val = sum(median_val)) %>%
    arrange(sum_median_val) %>%
    filter(parent_pesticide %in% unique(plot_dat$parent_pesticide))
  
  plot_dat$parent_pesticide <- factor(plot_dat$parent_pesticide, levels = order_chems$parent_pesticide)
  plot_dat$priority = TRUE
  # create full matrix of possibilities
  plot_dat_full <- expand.grid(unique(plot_dat$parent_pesticide), 
                               unique(plot_dat$metric_type)) %>%
    rename(parent_pesticide = Var1, metric_type = Var2) %>%
    left_join(plot_dat) %>%
    mutate(priority = ifelse(is.na(priority), FALSE, TRUE))
    
  
  p <- ggplot(plot_dat, aes(y = parent_pesticide, x = metric_type)) +
    geom_tile(fill = NA) +
    geom_text(aes(label = ifelse(priority, '*', '')), color = 'black', size = 7) +
    # scale_fill_viridis(option = 'viridis', direction = -1, begin = 0.1,
    #                    guide = guide_legend(direction = 'horizontal', 
    #                                         title.position = 'top', 
    #                                         legend.position = 'bottom', 
    #                                         label.position = 'bottom')) +
    facet_grid(rows = 'Class', scales = 'free', space = 'free_y') +
    scale_x_discrete(breaks = as.character(levels(plot_dat$metric_type)), labels = 
                       c(expression(paste("EAR"["mix"])),
                         expression(paste("EAR"["chem"])),
                         expression(paste("TQ"["chem"])),
                         expression(paste("Occurrence")))) +  
    labs(y = '', x = '') +
    theme_bw()+
    theme(legend.position = 'bottom',
          legend.spacing.x = unit(0.06, 'cm'), 
          strip.background = element_blank(),
          panel.border = element_rect(color = 'black'),
          panel.grid = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(0,0,0,0))
  
  #ggsave(filename = out_file, plot = p, height = 5, width = 2.8)
  chem_name_order <- ungroup(conc_median) %>%
    filter(parent_pesticide %in% keep_chems$parent_pesticide) %>%
    arrange(parent_pesticide, median_val) %>% mutate(order = 1:n()) %>%
    select(CAS, order)
    
  # which values were estimated
  bench_true <- raw_bench %>%
    select(chnm) %>%
    distinct()
  
  bench_max <- bench %>%
    group_by(chnm, CAS, parent_pesticide, site, date) %>%
    summarize(maxval = max(EAR)) %>% ungroup()
  
  all_samples <- bench %>%
    select(parent_pesticide, site, date) %>%
    distinct()
  
  all_compounds <- bench %>% 
    select(CAS, chnm, parent_pesticide) %>%
    distinct()
  
  all_possible <- full_join(all_samples, all_compounds)
  
  bench_median <- left_join(all_possible, ungroup(bench_max)) %>%
    mutate(maxval = ifelse(is.na(maxval), 0, maxval)) %>%
    group_by(CAS, chnm, parent_pesticide) %>%
    summarize(median_val = median(maxval)) %>%
    mutate(metric = 'TQchem',
           estimated = ifelse(chnm %in% bench_true$chnm, FALSE, TRUE)) %>%
    ungroup()

  plot_dat2 <- bind_rows(conc_median, ear_median, bench_median) %>%
    filter(parent_pesticide %in% keep_chems$parent_pesticide) %>%
    left_join(meta) %>% rename(parent_class = Class) %>%
    left_join(select(metadata, CAS, Class)) %>%
    mutate(degradate = grepl('Deg', Class)) %>%
    left_join(chem_name_order) %>%
    arrange(metric, parent_pesticide, order) %>% ungroup()
  
  plot_dat2$parent_class <- factor(plot_dat2$parent_class, levels = c('Herbicide', 'Fungicide', 'Insecticide'))
  plot_dat2$metric <- factor(plot_dat2$metric, levels = c('EARchem', 'TQchem', 'Concentration'),
                             labels = c(expression(paste("EAR"["chem"])), 
                                        expression(paste("TQ"["chem"])),
                                        expression(paste("Concentration [ug/L]"))))
  plot_dat2$parent_pesticide <- factor(plot_dat2$parent_pesticide, 
                                       levels = order_chems$parent_pesticide, labels = order_chems$parent_pesticide)
  
  #plot_dat2 <- mutate(plot_dat2, n = ifelse(degradate, NA, n))
  #plot_dat2$chnm <- factor(plot_dat2$chnm, levels = as.character(chem_name_order))
  
  #plot_dat2 <- arrange(plot_dat2, metric, parent_pesticide, )
  plot_dat2$n[plot_dat2$parent_pesticide %in% 'Chlorothalonil' & plot_dat2$metric %in% 'paste(\"Concentration [ug/L]\")'] <- 31
  text_dat <- group_by(plot_dat2, parent_pesticide, parent_class, metric) %>%
    summarize(median_val = sum(median_val) + ifelse(grepl('Conc', unique(metric)), 0.03, 0),
              n = mean(n, na.rm = TRUE))
  
  line_dat <- data.frame(metric = levels(plot_dat2$metric), median_val = c(0.001, 0.1, NA))

  my_label_parsed <- function (variable, value) {
    if (variable == "parent_pesticide") {
      return(as.character(value))
    } else {
      plyr::llply(as.character(value), function(x) parse(text = x))    
    }
  }

  p2 <- ggplot(data=plot_dat2, aes(x=parent_pesticide, y=median_val)) + 
    geom_bar(aes(color = estimated, fill = degradate), width=.8, size=.3, stat='identity') +  
    geom_text(data = text_dat, aes(label = n), color = 'red', size = 3, position = position_stack()) +
    coord_flip() +
    facet_grid(parent_class~metric, 
               space="free_y", scales="free", switch = 'x',
               labeller = my_label_parsed) +
    geom_hline(data = line_dat, aes(yintercept = median_val), color = 'blue') +
    theme_bw() +
    scale_color_manual(values = c('black', 'red'), name = 'Estimated', guide = guide_legend(title.position = 'top', label.position = 'bottom',
                                                                                                override.aes = list(fill = 'white'))) +
    scale_fill_manual(values = c('gray40', 'gray80'), name = 'Degradate', guide = guide_legend(title.position = 'top', label.position = 'bottom')) +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=.5),
          axis.text.x = element_text(size=8),
          legend.position = 'bottom', 
          legend.direction = 'horizontal',
          axis.text.y = element_blank(),
          legend.background = element_rect(fill = 'white', color = NA),
          strip.background = element_blank(),
          strip.placement = 'outside',
          legend.key.width = unit(1, 'cm'), plot.margin = margin(0,0,0,0)) +
    labs(x=NULL, y='', fill = 'Degradate', color = 'Estimated')
  
  p3 <- cowplot::get_legend(p + theme(legend.box.margin = margin(0,0,0,0)))
  p4 <- cowplot::get_legend(p2 + theme(legend.box.margin = margin(0,0,0,0)))
  
  
  row_1 <- plot_grid(p + theme(legend.position="none", plot.margin = margin(0,0,0,0)),
            p2 + theme(legend.position= 'none', plot.margin = margin(0,0,0,0)),
            nrow = 1, rel_widths = c(1,3), align = 'h', axis = 'bt')
  row_2 <- plot_grid(p3, p4, nrow = 1, rel_widths = c(2.2,3))
 
  cowplot::plot_grid(p3, p4, row_1, nrow = 2, rel_heights = c(0.1, 0.1, 1), rel_widths = c(0.1, 0.1, 10))
  
  ggsave(filename = out_file, 
         plot_grid(row_1, row_2,rel_heights = c(10,1.5),rel_widths = c(1, 0.5), align = 'h', ncol = 1),
         height = 5, width = 8)
    
  
  
  }
  



  

