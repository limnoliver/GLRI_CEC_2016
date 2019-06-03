# plot DL values
# these are median values of DLs
# order chems by ...

plot_dls_tox <- function(chem_dls, chem_dl_tox, pcodeinfo, fig_name){
  
  chem_info <- select(chem_dl_tox, CAS, chnm) %>%
    distinct() %>%
    left_join(filter(select(pcodeinfo, CAS = casrn, pCode = parameter_cd), pCode %in% chem_dls$pCode))
  
  #chem_info <- make('chem_info')
  
  plot_data_dl <- toxEval::tox_boxplot_data(chem_dl_tox, category = 'Chemical', sum_logic = TRUE)
  
  plot_data_dl <- left_join(plot_data_dl, chem_info) %>%
    left_join(chem_dls) %>%
    mutate(prop_detected = n_detect/n_all)

  
  p <- ggplot(plot_data_dl, aes(x = meanEAR, y = chnm)) +
    geom_point(aes(color = Class, size = prop_detected), alpha = 0.8) +
    scale_x_log10() +
    geom_vline(xintercept = 0.001) +
    labs(size = 'Proportion of samples \nwith detections', y = '') +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7))
  
  ggsave(fig_name, height = 10, width = 7)
}

