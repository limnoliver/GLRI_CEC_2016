# ranking figs

make_site_fig <- function(file_name, site_rankings, sites) {
  
  site_dat <- select(sites, site = site_no, landuse = `Dominant.land.use.`, shortName, 
                     ag_total = `Ag..total`, Urban) %>%
    mutate(landuse = ifelse(landuse == 'Wetland' | landuse == 'Forest', "Natural", landuse),
           perc_ag_urban = ag_total + Urban)
  
  plot_dat <- site_rankings
  plot_dat$variable <- factor(plot_dat$variable,
                              levels = c('n_detected', 'n_detected_new',
                                         'mean_conc', 'sum_conc',
                                         'mean_bench', 'max_bench',
                                         'mean_sumEAR', 'max_sumEAR'))
  var_summary <- group_by(plot_dat, site, variable) %>%
    summarize(y0 = quantile(value, 0.1),
              y25 = quantile(value, 0.25),
              y50 = median(value),
              y75 = quantile(value, 0.75), 
              y100 = quantile(value, 0.9), 
              rank = mean(max_rank)) %>%
    left_join(site_dat)
  
  var_summary$landuse <- factor(var_summary$landuse, levels = c('Natural','Urban', 'AgMix', 'Crops'))
  var_summary <- arrange(var_summary, landuse, perc_ag_urban)
  site_order <- unique(var_summary$shortName)
  var_summary$shortName <- factor(var_summary$shortName, levels = site_order)
  
  var_ymax <- var_summary %>%
    group_by(variable) %>%
    summarize(ylim = max(y50)*1.2)
  
  variable_dat <- filter(plot_dat, site != '04193500') %>%
    group_by(site, variable) %>%
    summarize(ymax_quart = (quantile(value, .75) * 1.7), 
              ymax_obs = max(value)) %>%
    group_by(variable) %>%
    summarize(ymax_quart = max(ymax_quart), 
              ymax_obs = max(ymax_obs)) %>%
    rowwise() %>%
    mutate(ymax = min(ymax_quart, ymax_obs))
  
  # plot_dat <- filter(plot_dat, site != '04193500')
  vars <- unique(plot_dat$variable)
  var_ylab <- c("No. chemicals detected", "No. chemicals unique to sample", "Sum Concentration", 
                "Mean concentration", "Mean of sumEAR", "Max of sumEAR", 
                "Mean benchmark ratio", "Max benchmark ratio")
  p <- list()
  q <- list()
  for (i in 1:length(vars)) {
    temp_dat <- filter(plot_dat, variable == vars[i])
    p[[i]] <- ggplot(temp_dat, aes(x = site, y = value)) +
      geom_boxplot(aes(fill = max_rank)) +
      #facet_wrap(~variable, nrow = 4, scales = 'free_y') +
      scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                           midpoint = 8.5, limits = c(1, 16)) +
      coord_cartesian(ylim = c(0, variable_dat$ymax_quart[variable_dat$variable == vars[i]])) +
      labs(x = '', y = var_ylab[i], fill = "Site Rank") +
      theme(axis.text.x = element_blank())
    
    temp_dat <- filter(var_summary, variable == vars[[i]])
    q[[i]] <- ggplot(temp_dat, aes(x = shortName)) +
      geom_boxplot(aes(fill = rank, ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                   stat = 'identity') +
      #facet_wrap(~variable, nrow = 4, scales = 'free_y') +
      scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                           midpoint = 8.5, limits = c(1, 16), breaks = c(1,4,8,12,16)) +
      #coord_cartesian(ylim = c(0, var_ymax$ylim[var_ymax$variable == vars[i]])) +
      labs(x = '', y = '', title = var_ylab[i], fill = "Site Rank") +
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.05, vjust = 0.1, size = 16),
            plot.margin = margin(0, 2, 0, 0))
    
    if (i == 6)
     p[[i]] <- p[[i]] +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = 'Site')
  }
  
  # modify four figures to use for cowplot fig
  q1 <- q[[1]] + theme(legend.position = c(0.4, .9),
                       legend.direction = 'horizontal', 
                       legend.key.height = unit(0.5, 'cm'),
                       legend.margin = margin(5,7,5,5),
                       legend.box.background = element_rect(color = 'black'))
  q2 <- q[[3]] + theme(legend.position = 'none') +
    scale_y_continuous(trans = 'log10')
  q3 <- q[[8]] + theme(legend.position = 'none') +
    scale_y_continuous(trans = 'log10')
  q4 <- q[[6]] + theme(legend.position = 'none',
                       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_y_continuous(trans = 'log10')
  
  legend_plot <- q4 + theme(legend.position = 'bottom', 
                            legend.key.height = unit(0.5, 'cm'),
                            legend.box.margin = margin(-25,0,0,0), 
                            legend.margin = margin(0,0,0,0), 
                            legend.direction = 'horizontal')
  
  # if you want to go back to old figs, just uncomment the coord_cartesian line
  # and get rid of transformation in plot post-processing.
  # think if we use raw data instead of summary statistics, would look better in 
  # log space.
  legend_b <- get_legend(legend_plot)
  qrow <- plot_grid(q1, q2, q3, q4, ncol = 1, rel_heights = c(1,1,1,1.3),
                    align = 'v', axis = 'l')
  
  ggsave(file_name, qrow, height = 12, width = 9)
}

make_site_tile <- function(file_name, site_rankings, sites) {
  
  site_dat <- select(sites, site = site_no, landuse = `Dominant.land.use.`, shortName, 
                     ag_total = `Ag..total`, ag_crops = `Ag..crops`, Urban) %>%
    mutate(landuse = ifelse(landuse == 'Wetland' | landuse == 'Forest', "Natural", landuse),
           perc_ag_urban = ag_total + Urban,
           perc_row_urban = ag_crops+Urban)
  
  plot_dat <- site_rankings
  plot_dat$metric <- factor(plot_dat$metric)
  levels(plot_dat$metric) <- c(levels(plot_dat$metric)[1], 'med hits/sample',levels(plot_dat$metric)[3], 
                               'max TQchem', 'max chems detected', 'med hits/sample', 
                               'med TQchem', 'med detects/sample', 'months w/hits',
                               'total chems detected', 'months w/hits')
  site_dat$landuse <- factor(site_dat$landuse, levels = rev(c('Natural','Urban', 'AgMix', 'Crops')))
  site_dat <- arrange(site_dat, landuse, perc_ag_urban)
  site_order <- unique(site_dat$shortName)
  site_dat$shortName <- factor(site_dat$shortName, levels = site_order)
  
  # filter out specific vars that we don't want to use
  #plot_dat <- filter(plot_dat, variable %in% c('No. chems detected', 'Mean concentration',
  #                                             'Mean benchmark ratio', 'Mean sumEAR'))
  
  plot_dat <- left_join(plot_dat, site_dat)
  
  
  plot_dat_min <- plot_dat %>%
    group_by(metric, metric_type) %>%
    summarize(metric_value = metric_value[which.min(metric_value)]) %>%
    mutate(print_min = metric_value)
  
  plot_dat <- plot_dat %>%
    mutate(print_max = ifelse(relative_value == 1, metric_value, NA)) 
  plot_dat <- left_join(plot_dat, plot_dat_min, by = c("metric", "metric_value", "metric_type", 'metric_value'))%>%
    mutate(print_min = ifelse(is.na(print_min), '', format(print_min, scientific = FALSE, digits = 4)))
  
  
  #long_plot_dat <- gather(plot_dat, key = rank_variable, value = rank_value, median_rank, max_rank) %>%
  #  mutate(rank_variable = ifelse(rank_variable == 'max_rank', "Acute (max) rank", "Chronic (med) rank"))
  
  p <- ggplot(plot_dat, aes(y = shortName, x = metric)) +
    geom_tile(aes(fill = relative_value), color = 'black') +
    geom_text(aes(label = round(print_max, 1)), color = 'white') +
    geom_text(aes(label = print_min), color = 'black') +
    scale_fill_viridis(option = 'viridis', direction = -1, begin = 0.1) +
    #facet_wrap(landuse~rank_variable, ncol = 2, scales = 'free_y') +
    facet_grid(landuse~metric_type, scales = 'free', space = 'free_y') +
    labs(y = '', x = '', fill = "Relative Value") +
    theme(axis.text.x = element_text(angle = 45, v = 1, h = 1),
          legend.position = 'bottom', legend.direction = 'horizontal',
          legend.title = element_text(vjust = 1.5, hjust = 1.5))
  
  ggsave(file_name, p, height = 7, width = 10)
}
add_final_rank <- function(site_rankings, sites) {
  # create wide version of site rank data
  wide_dat <- site_rankings
  wide_dat$metric <- factor(wide_dat$metric)
  levels(wide_dat$metric) <- c(levels(wide_dat$metric)[1], 'med hits/sample',levels(wide_dat$metric)[3], 
                               'max TQchem', 'max chems detected', 'med hits/sample', 
                               'med TQchem', 'med detects/sample', 'months w/hits',
                               'total chems detected', 'months w/hits')
  
  wide_dat_values <- wide_dat %>%
    ungroup() %>%
    mutate(column_name = paste(metric_type, metric, sep = '_')) %>%
    select(-metric, -rank_value, -relative_value, -metric_type) %>%
    tidyr::spread(column_name, metric_value)
  
  wide_rel_vals <- wide_dat %>%
    ungroup() %>%
    mutate(column_name = paste(metric_type, metric, sep = '_')) %>%
    select(-metric, -metric_value, -rank_value, -metric_type) %>%
    tidyr::spread(column_name, relative_value)
  
  # calculate mean values -- and weight benchmarks by 3 (most confident in those #s), 
  # EARs by 2 (better than just a chem being present, but not exactly sure how it scales), 
  # and occurence data by 1
  rel_vals <- wide_rel_vals %>%
    mutate(bench_mean_rel = round(rowMeans(select(., starts_with('Bench')), na.rm = TRUE),2),
           earmix_mean_rel = round(rowMeans(select(., starts_with('EARmix')), na.rm = TRUE),2),
           occurance_mean_rel = round(rowMeans(select(., starts_with('Occurance')), na.rm = TRUE),2))
  
  rel_vals$final_rel_val_321 <- round((3*rel_vals$bench_mean_rel + 2*rel_vals$earmix_mean_rel + rel_vals$occurance_mean_rel)/6, 2)
  #rel_vals$final_rel_val_111 <- (rel_vals$bench_mean_rel + rel_vals$earmix_mean_rel + rel_vals$occurance_mean_rel)/3

  dat_add <- rel_vals %>%
    select(site, bench_mean_rel, earmix_mean_rel, occurance_mean_rel, final_rel_val_321) %>%
    tidyr::gather(key = 'metric', value = 'relative_value', -site) %>%
    mutate(metric_type = case_when(
      grepl('bench', metric) ~ 'Bench',
      grepl('earmix', metric) ~ 'EARmix',
      grepl('occurance', metric) ~ 'Occurance',
      grepl('final', metric) ~ 'RRI'
    ))
  
  dat_out <- bind_rows(wide_dat, dat_add)
  
  return(dat_out)
}
create_site_rank_table <- function(site_rankings_final, sites, out_file) {
  
 rel_vals <- filter(site_rankings_final, grepl('mean_rel|final_rel', metric)) %>%
   ungroup(metric) %>%
   select(-metric_type, -rank_value, -metric_value) %>%
   tidyr::spread(key = metric, value = relative_value)
 
 values <- filter(site_rankings_final, !grepl('mean_rel|final_rel', metric)) %>%
   ungroup(metric) %>%
   mutate(metric = paste(metric_type, metric, sep = '_')) %>%
   select(-metric_type, -rank_value, -relative_value) %>%
   tidyr::spread(key = metric, value = metric_value)
 
 dat_out <- left_join(rel_vals, values) %>%
   left_join(select(sites, site, shortName)) %>%
   select(site, starts_with('Bench', ignore.case = TRUE), 
          starts_with('EAR', ignore.case = TRUE), 
          starts_with('Occurance', ignore.case = TRUE), 
          RRI = final_rel_val_321) %>%
   arrange(-RRI) %>%
   
  
  write.csv(dat_out, out_file, row.names = FALSE)
  
}

plot_rankings_v_landuse <- function(){
  plot_vals <- select(bench_rel_vals, site, final_rel_val_321, final_rel_val_111) %>%
    left_join(site_dat)
  
  ggplot(plot_vals, aes(y=final_rel_val_111, x = perc_row_urban)) +
    geom_point(aes(color = landuse))
}

