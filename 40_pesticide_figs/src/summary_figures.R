sum_pest_conc <- function(reduced_dat, chems_missing_cas) {
  dat <- select(reduced_dat, SiteID, pCode, value, remark_cd, sample_dt, source, State, Site, Glyphosate, Neonics) %>%
    filter(!(pCode %in% chems_missing_cas)) %>%
    distinct()
  sum_chems <- dat %>%
    #mutate(date_no_time = date(pdate)) %>%
    group_by(SiteID, sample_dt) %>%
    summarize(sum_conc = sum(value, na.rm = T), 
              sum_conc_detect = sum(value[which(!(remark_cd %in% "<"))], na.rm = T),
              n_chems = n(), 
              neonic_meas = ifelse('neonic' %in% source, TRUE, FALSE),
              n_detected = length(which(!(remark_cd %in% "<"))))
  
  sum_detected_chems <- filter(dat, !(remark_cd %in% "<")) %>%
    group_by(pCode) %>%
    summarise(n_detected = n(),
              mean_val = mean(value),
              median_val = median(value),
              min_val = min(value))

  sum_study <- ungroup(sum_chems) %>%
    summarise_at(vars(sum_conc_detect, n_detected), funs(mean, median, min, max, sd))
  
  
  sum_sites <- sum_chems %>%
    group_by(SiteID) %>%
    summarise_at(vars(sum_conc_detect, n_detected), funs(mean, median, min, max, sd))
  
  #test <- mutate(merged_dat, date_no_time = date(pdate)) %>%
  # filter(SiteID %in% '04193500' & date_no_time %in% as.Date('2016-06-07'))
  
  
  # basicaly get rid of single sample with n = 1 chem
  # need to fix at some point
  #sum_chems_subset <- filter(sum_chems, n_chems >220)
  
   sum_chems_subset_order <- group_by(sum_chems, SiteID) %>%
     summarize(median = median(sum_conc)) %>%
     arrange(median) %>%
     pull(SiteID)
  
  sum_chems_subset_order_detect <- group_by(sum_chems, SiteID) %>%
    summarize(median = median(sum_conc_detect)) %>%
    arrange(median) %>%
    pull(SiteID)
  
  # miscellaneous summaries used in past to detect issues in data:
  # names_chems = paste(chnm, collapse = ', ')
  # n_missing_chems = length(which(!(all_chems %in% chnm)))
  # names_missing_chems = paste(all_chems[which(!(all_chems %in% chnm))], collapse = ', ')
  
  sum_chems$SiteID_all <- factor(sum_chems$SiteID, levels = sum_chems_subset_order)
  sum_chems$SiteID_detected <- factor(sum_chems$SiteID, levels = sum_chems_subset_order_detect)
  # 
  # 
  # site_summ <- group_by(sum_chems, SiteID) %>%
  #   summarize(n = n())
  
  sum_chems$month <- month(sum_chems$sample_dt, label = TRUE)
  sum_chems$sum_conc_ppb <- sum_chems$sum_conc/1000
  sum_chems$sum_conc_detect_ppb <- sum_chems$sum_conc_detect/1000
  
  
  return(sum_chems)
}

combine_dat <- function(sum_dat, site_dat) {
  graph_dat_thresh <- make('parent_sums') %>%
    left_join(make('parent_class')) %>%
    mutate(day = lubridate::yday(date)) %>%
    mutate(month = lubridate::month(date)) %>%
    filter(type == 'p_sumval' & measure_type %in% c('ear', 'bench')) %>%
    filter(!is.na(sumval)) %>%
    group_by(site, month, measure_type) %>%
    summarize(max_val = max(sumval, na.rm = TRUE)) %>%
    ungroup()
  
  graph_dat_hits <- make('parent_sums') %>%
    left_join(make('parent_class')) %>%
    mutate(day = lubridate::yday(date)) %>%
    mutate(month = lubridate::month(date)) %>%
    filter(type == 'p_sumval'  & measure_type %in% c('ear', 'bench')) %>%
    filter(measure_type == 'ear' & sumval >= 0.001 | measure_type == 'bench' & sumval >= 0.01) %>%
    group_by(site, month, day, measure_type) %>%
    summarize(n_hits = n()) %>%
    ungroup()
  
  maumee_fix <- filter(graph_dat_hits, site == '04193500') %>%
    group_by(site, month, measure_type) %>%
    summarize(n_hits = max(n_hits))
    
  graph_dat_hits <- filter(graph_dat_hits, !site == '04193500') %>%
    bind_rows(maumee_fix)
    
  # 
  # med_graph_dat <- graph_dat_thresh %>%
  #   group_by(month, measure_type) %>%
  #   summarize(median_max_val = median(max_val))
  
  full_matrix <- expand.grid(site = unique(graph_dat_hits$site), month = 1:12, measure_type = c('bench', 'ear'))
  
  hits <- left_join(full_matrix, graph_dat_hits) %>%
    mutate(n_hits = ifelse(is.na(n_hits), 0, n_hits))
  
  thresh <- left_join(full_matrix, graph_dat_thresh) %>%
    mutate(max_val = ifelse(is.na(max_val), min(graph_dat_thresh$max_val)*0.5, max_val))
  
  # med_graph_dat_hit <- hits %>%
  #   group_by(month, measure_type) %>%
  #   summarize(n_hits = median(n_hits), site = 'all')
  # 
  # all_dat_hits <- bind_rows(hits, med_graph_dat_hit) %>%
  #   mutate(site_cat = ifelse(site %in% 'all', 'all', 'individual'))
  
  
  sites <- make('sites') %>%
    select(site = site_no, dominant_lu = Dominant.land.use.) %>%
    mutate(dominant_lu = ifelse(dominant_lu %in% c('Wetland', 'Forest'), 'Wetland/Forest', dominant_lu))
    
  all_dat <- left_join(graph_dat, sites) %>%
    add_row(month = med_graph_dat$month, max_val = med_graph_dat$median_max_val, site = 'all', measure_type = med_graph_dat$measure_type) %>%
    mutate(site_cat = ifelse(site %in% 'all', 'all', 'individual'))
  
  library(ggplot2)
  line_dat <- data.frame(x = 1, y = 0.01, measure_type = c('bench', 'ear'), z = c(0.1, 0.001))
  ggplot(data = thresh, aes(x = month, y = max_val)) +
    geom_boxplot(aes(group = month)) +
    #geom_point() +
    #geom_line(aes(group = site, color = site_cat), alpha = 0.5) +
    #scale_color_manual(values = c('black', 'gray')) +
    scale_y_log10() +
    facet_wrap(~measure_type) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  ggplot(data = hits, aes(x = month, y = n_hits)) +
    geom_boxplot(aes(group = month)) +
    #geom_point(aes( color = site_cat), alpha = 0.5) +
    #geom_line(aes(group = site, color = site_cat), alpha = 0.5) +
    #scale_color_manual(values = c('black', 'gray')) +
    #scale_y_log10() +
    facet_wrap(~measure_type) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  }


boxplot_bysite <- function(sum_conc, target_name, detect_only = TRUE) {
  if (detect_only == FALSE) {
    p <- ggplot(sum_conc, aes(x = SiteID_all, y = sum_conc_ppb)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
      geom_jitter(aes(shape=neonic_meas), position=position_jitter(0.1), alpha = 0.9, col = 'darkgray') +
      scale_shape_manual(values = c(1, 16), name = 'Neonics measured') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, h = 1)) +
      labs(y = 'Total pesticide concentration (ppb)')
    
  } else {
    p <- ggplot(sum_conc, aes(x = SiteID_detected, y = sum_conc_detect_ppb)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
      geom_jitter(aes(shape=neonic_meas), position=position_jitter(0.1), alpha = 0.9, col = 'darkgray') +
      scale_shape_manual(values = c(1, 16), name = 'Neonics measured') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, h = 1)) +
      labs(y = 'Total pesticide concentration (ppb)') +
      scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 5, 10, 20, 30, 40, 50))
  }
  
  # if (by_chem_group == TRUE) {
  #   p <- p +
  #     facet_wrap()
  # }
  
  
  ggsave(target_name, p, height = 4, width = 6)
  
  
}

boxplot_bysite_month <- function(sum_conc, target_name, site_dat) {
  
  dat <- left_join(sum_conc, site_dat, by = 'SiteID')
  
  ordered_names <- dat %>%
    ungroup() %>%
    select(`Short Name`, SiteID_detected) %>%
    arrange(SiteID_detected) %>%
    distinct()
  
  dat$`Short Name` <- factor(dat$`Short Name`, levels = ordered_names$`Short Name`)

  
  month.cols<-c(viridis(6, begin=.2, end=.99), rev(magma(6, begin=.2, end=.95)))
  p_month <- ggplot(dat, aes(x = `Short Name`, y = sum_conc_detect_ppb)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_jitter(aes(shape=neonic_meas, color = month), position=position_jitter(0.15), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = month.cols, name = 'Month') +
    scale_shape_manual(values = c(17, 16), name = 'Neonics measured') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, h = 1)) +
    labs(y = 'Total pesticide concentration (ppb)', x = 'Site') +
    scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 5, seq(10, 50, 10))) +
    guides(color = guide_legend(ncol = 2),shape = guide_legend(ncol = 2))
  
  ggsave(target_name, p_month, height = 4, width = 6)
  
  
}

boxplot_ndetect_bysite_month <- function(sum_conc, target_name, site_dat) {
  dat <- left_join(sum_conc, site_dat, by = 'SiteID')
  
  ordered_names <- dat %>%
    ungroup() %>%
    select(`Short Name`, SiteID_detected) %>%
    arrange(SiteID_detected) %>%
    distinct()
  
  dat$`Short Name` <- factor(dat$`Short Name`, levels = ordered_names$`Short Name`)
  
  month.cols<-c(viridis(6, begin=.2, end=.99), rev(magma(6, begin=.2, end=.95)))
  p_month <- ggplot(dat, aes(x = `Short Name`, y = n_detected)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_jitter(aes(shape=neonic_meas, color = month), position=position_jitter(0.15), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = month.cols, name = 'Month') +
    scale_shape_manual(values = c(17, 16), name = 'Neonics measured') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, h = 1)) +
    labs(y = 'Number of pesticides detected', x = 'Sites from low to high median concentration') +
    #scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 5, seq(10, 50, 10))) +
    guides(color = guide_legend(ncol = 2),shape = guide_legend(ncol = 2))
  
  ggsave(target_name, p_month, height = 4, width = 6)
  
  
}

boxplot_bymonth_site <- function(sum_conc, target_name, site_dat) {
  
  dat <- left_join(sum_conc, site_dat, by = 'SiteID')
  
  ordered_names <- dat %>%
    ungroup() %>%
    select(`Short Name`, SiteID_detected) %>%
    arrange(SiteID_detected) %>%
    distinct()
  
  dat$`Short Name` <- factor(dat$`Short Name`, levels = ordered_names$`Short Name`)
  
  site.cols<- inferno(16, begin=0.4, end=.9, direction = -1)
  p <- ggplot(dat, aes(x = month, y = sum_conc_detect_ppb)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_jitter(aes(color = `Short Name`, size = n_detected), position=position_jitter(0.15), alpha = 0.4) +
    scale_color_manual(values = site.cols, name = 'Sites (ranked low to high \nmedian concentration)') +
    #scale_shape_manual(values = c(17, 16), name = 'Neonics measured') +
    scale_size_continuous(range = c(0.7, 5), breaks = c(0, 5, 10, 20, 40, 60), name = '# chemicals detected') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, h = 1)) +
    labs(y = 'Total concentration of \ndetected pesticides (ppb)', x = '') +
    guides(color = guide_legend(ncol = 2),size = guide_legend(ncol = 2)) +
    scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 5, seq(10, 50, 10)))
  
  ggsave(target_name, p, height = 4, width = 6)
  
  
}

boxplot_ndetect_bymonth_site <- function(sum_conc, target_name, site_dat) {
  dat <- left_join(sum_conc, site_dat, by = 'SiteID')
  
  ordered_names <- dat %>%
    ungroup() %>%
    select(`Short Name`, SiteID_detected) %>%
    arrange(SiteID_detected) %>%
    distinct()
  
  dat$`Short Name` <- factor(dat$`Short Name`, levels = ordered_names$`Short Name`)
  
  site.cols<- inferno(16, begin=0.2, end=.9, direction = -1)
  p <- ggplot(dat, aes(x = month, y =n_detected)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_jitter(aes(color = `Short Name`, size = sum_conc_detect_ppb), position=position_jitter(0.15), alpha = 0.4) +
    scale_color_manual(values = site.cols, name = 'Sites (ranked low to high \nmedian concentration)') +
    #scale_shape_manual(values = c(17, 16), name = 'Neonics measured') +
    scale_size_continuous(range = c(0.7, 7), breaks = c(0, 1, 5, 10, 20, 40), name = 'Concentration of detected \npesticides (ppb)') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, h = 1)) +
    labs(y = 'Number of pesticides detected', x = '') +
    guides(color = guide_legend(ncol = 2),size = guide_legend(ncol = 2))
  
  ggsave(target_name, p, height = 5, width = 7)
  
  
}

plot_throughtime_bysite <- function(sum_conc, target_name, site_dat) {
  
  dat <- left_join(sum_conc, site_dat, by = 'SiteID')
  
  ordered_names <- dat %>%
    ungroup() %>%
    select(`Short Name`, SiteID_detected) %>%
    arrange(SiteID_detected) %>%
    distinct()
  
  dat$`Short Name` <- factor(dat$`Short Name`, levels = ordered_names$`Short Name`)
  
  p2 <- ggplot(dat, aes(x = sample_dt, y = sum_conc_detect_ppb)) +
    geom_point(aes(size = n_detected), alpha = 0.5) +
    scale_size_continuous(range = c(0.7, 5), breaks = c(0, 5, 10, 20, 40, 60), name = '# chemicals detected') +
    facet_wrap(~`Short Name`, ncol = 2, scales = 'free_y') +
    theme_bw() +
    labs(x = '', y = 'Total pesticide concentration (ppb)') +
    scale_y_continuous(expand = expand_scale(mult = c(.15, .25))) +
    theme(axis.text.x = element_text(angle = 45, h = 1))
  
  ggsave(target_name, p2, height = 8, width = 6)
  
}

plot_nchem_throughtime_bysite <- function(sum_conc, target_name, site_dat) {
  
  dat <- left_join(sum_conc, site_dat, by = 'SiteID')
  
  ordered_names <- dat %>%
    ungroup() %>%
    select(`Short Name`, SiteID_detected) %>%
    arrange(SiteID_detected) %>%
    distinct()
  
  dat$`Short Name` <- factor(dat$`Short Name`, levels = ordered_names$`Short Name`)
  
  p2 <- ggplot(dat, aes(x = sample_dt, y = n_detected)) +
    geom_point(aes(size = sum_conc_detect_ppb), alpha = 0.5) +
    scale_size_continuous(range = c(1, 6), breaks = c(0, 5, 10, 20, 40), name = 'Concentration of detected \nchemicals (ppb)') +
    facet_wrap(~`Short Name`, ncol = 2, scales = 'free_y') +
    theme_bw() +
    labs(x = '', y = 'Number of pesticides detected') +
    scale_y_continuous(expand = expand_scale(mult = c(.15, .25))) +
    theme(axis.text.x = element_text(angle = 45, h = 1))
  
  ggsave(target_name, p2, height = 8, width = 8)
  
}

plot_nchem_bysitedate <- function(reduced_dat, target_name) {
  
  count_dat <- group_by(reduced_dat, SiteID, pCode) %>%
    summarise(count = n())
  
  missing_pcodes <- spread(count_dat, key = SiteID, value = count) %>%
    rowwise() %>%
    mutate(n_na = sum(is.na()))
  
  ggplot(count_dat, aes(x = SiteID, y = pCode)) +
    geom_tile(aes(fill = factor(count))) 
    
}
