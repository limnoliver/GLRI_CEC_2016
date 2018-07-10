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
  
  # sum_detected_chems <- filter(dat, !(remark_cd %in% "<")) %>%
  #   group_by(pCode) %>%
  #   summarise(n_detected = n(), 
  #             mean_val = mean(value), 
  #             median_val = median(value), 
  #             min_val = min(value))
  
  # sum_study <- ungroup(sum_chems) %>%
  #   summarise_at(vars(sum_conc_detect, n_detected), funs(mean, median, min, max, sd))
  # 
  #   
  # sum_sites <- sum_chems %>%
  #   group_by(SiteID) %>%
  #   summarise_at(vars(sum_conc_detect, n_detected), funs(mean, median, min, max, sd))
  # 
  # #test <- mutate(merged_dat, date_no_time = date(pdate)) %>%
  # #  filter(SiteID %in% '04193500' & date_no_time %in% as.Date('2016-06-07'))
  # 
  # 
  # # basicaly get rid of single sample with n = 1 chem
  # # need to fix at some point
  # #sum_chems_subset <- filter(sum_chems, n_chems >220)
  # sum_chems_subset_order <- group_by(sum_chems_subset, SiteID) %>%
  #   summarize(median = median(sum_conc)) %>%
  #   arrange(median) %>%
  #   pull(SiteID)
  # 
  # sum_chems_subset_order_detect <- group_by(sum_chems_subset, SiteID) %>%
  #   summarize(median = median(sum_conc_detect)) %>%
  #   arrange(median) %>%
  #   pull(SiteID)
  # 
  # # miscellaneous summaries used in past to detect issues in data:
  # #names_chems = paste(chnm, collapse = ', '),
  # #n_missing_chems = length(which(!(all_chems %in% chnm))),
  # #names_missing_chems = paste(all_chems[which(!(all_chems %in% chnm))], collapse = ', ')  
  # sum_chems_subset$SiteID_all <- factor(sum_chems_subset$SiteID, levels = sum_chems_subset_order)
  # sum_chems_subset$SiteID_detected <- factor(sum_chems_subset$SiteID, levels = sum_chems_subset_order_detect)
  # 
  # 
  # site_summ <- group_by(sum_chems, SiteID) %>%
  #   summarize(n = n())
  
  sum_chems$month <- month(sum_chems$sample_dt, label = TRUE)
  sum_chems$sum_conc_ppb <- sum_chems$sum_conc/1000
  sum_chems$sum_conc_detect_ppb <- sum_chems$sum_conc_detect/1000
  
  
  return(sum_chems)
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
  
  if (by_chem_group == TRUE) {
    p <- p +
      facet_wrap()
  }
  
  
  ggsave(target_name, p, height = 4, width = 6)
  
  
}

boxplot_bysite_month <- function(sum_conc, target_name) {
  
  month.cols<-c(viridis(6, begin=.2, end=.99), rev(magma(6, begin=.2, end=.95)))
  p_month <- ggplot(sum_conc, aes(x = SiteID_detected, y = sum_conc_detect_ppb)) +
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

boxplot_ndetect_bysite_month <- function(sum_conc, target_name) {
  
  month.cols<-c(viridis(6, begin=.2, end=.99), rev(magma(6, begin=.2, end=.95)))
  p_month <- ggplot(sum_conc, aes(x = SiteID_detected, y = n_detected)) +
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

boxplot_bymonth_site <- function(sum_conc, target_name) {
  
  site.cols<- inferno(16, begin=0.4, end=.9, direction = -1)
  p <- ggplot(sum_conc, aes(x = month, y = sum_conc_detect_ppb)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_jitter(aes(color = SiteID_detected, size = n_detected), position=position_jitter(0.15), alpha = 0.4) +
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

boxplot_ndetect_bymonth_site <- function(sum_conc, target_name) {
  
  site.cols<- inferno(16, begin=0.2, end=.9, direction = -1)
  p <- ggplot(sum_conc, aes(x = month, y =n_detected)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_jitter(aes(color = SiteID_detected, size = sum_conc_detect_ppb), position=position_jitter(0.15), alpha = 0.4) +
    scale_color_manual(values = site.cols, name = 'Sites (ranked low to high \nmedian concentration)') +
    #scale_shape_manual(values = c(17, 16), name = 'Neonics measured') +
    scale_size_continuous(range = c(0.7, 7), breaks = c(0, 1, 5, 10, 20, 40), name = 'Concentration of detected \npesticides (ppb)') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, h = 1)) +
    labs(y = 'Number of pesticides detected', x = '') +
    guides(color = guide_legend(ncol = 2),size = guide_legend(ncol = 2))
  
  ggsave(target_name, p, height = 5, width = 7)
  
  
}

plot_throughtime_bysite <- function(sum_conc, target_name) {
  p2 <- ggplot(sum_conc, aes(x = sample_dt, y = sum_conc_detect_ppb)) +
    geom_point(aes(size = n_detected), alpha = 0.5) +
    scale_size_continuous(range = c(0.7, 5), breaks = c(0, 5, 10, 20, 40, 60), name = '# chemicals detected') +
    facet_wrap(~SiteID_detected, ncol = 2, scales = 'free_y') +
    theme_bw() +
    labs(x = '', y = 'Total pesticide concentration (ppb)') +
    scale_y_continuous(expand = expand_scale(mult = c(.15, .25))) +
    theme(axis.text.x = element_text(angle = 45, h = 1))
  
  ggsave(target_name, p2, height = 8, width = 6)
  
}

plot_nchem_throughtime_bysite <- function(sum_conc, target_name) {
  p2 <- ggplot(sum_conc, aes(x = sample_dt, y = n_detected)) +
    geom_point(aes(size = sum_conc_detect_ppb), alpha = 0.5) +
    scale_size_continuous(range = c(1, 6), breaks = c(0, 5, 10, 20, 40), name = 'Concentration of detected \nchemicals (ppb)') +
    facet_wrap(~SiteID_detected, ncol = 2, scales = 'free_y') +
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