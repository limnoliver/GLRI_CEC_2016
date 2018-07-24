# ranking figs

make_site_fig <- function(site_rankings) {

  plot_dat <- site_rankings
  plot_dat$variable <- factor(plot_dat$variable,
                              levels = c('n_detected', 'n_detected_new',
                                         'mean_conc', 'sum_conc',
                                         'mean_bench', 'max_bench',
                                         'mean_sumEAR', 'max_sumEAR'))
  
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
  p <- list()
  for (i in 1:length(vars)) {
    temp_dat <- filter(plot_dat, variable == vars[i])
    p[[i]] <- ggplot(temp_dat, aes (x = site, y = value)) +
      geom_boxplot(aes(fill = max_rank)) +
      #facet_wrap(~variable, nrow = 4, scales = 'free_y') +
      scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                           midpoint = 8) +
      coord_cartesian(ylim = c(0, variable_dat$ymax_obs[variable_dat$variable == vars[i]])) +
      labs(x = 'Site', y = vars[i], fill = "Site Rank")
  }
  # use precomputed statistics instead - see example in https://ggplot2.tidyverse.org/reference/geom_boxplot.html
  p <- ggplot(plot_dat, aes (x = site, y = value)) +
    geom_boxplot(aes(fill = max_rank)) +
    facet_wrap(~variable, nrow = 4, scales = 'free_y') +
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                        midpoint = 8)
  
  p_dat <- ggplot_build(p)

  for (i in 2:8) {
    p_dat$layout$panel_scales_y[[i]]$range$range <- c(0, variable_dat$ymax[i])
  }
  print(p_dat)
  test <- ggplot_gtable(p_dat)
  print(test)
  plot(test)
}