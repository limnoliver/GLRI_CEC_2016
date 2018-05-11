library(toxEval)

graphData_combos <- function(graphData_tox, graphData_wq, graphData_conc, file_to_save){
  

  graphData_tox$guide_side <- "ToxCast\nMaximum EAR per Site"
  graphData_wq$guide_side <- "Traditional\nMaximum Toxicity Quotient per Site"
  graphData_conc$guide_side <- "Concentration [ug/L]"

  toxPlot_wq <- combo_plot_matches(graphData_tox, graphData_wq, 
                                   thres_1 = NA, thres_2 = NA, 
                                   drop = FALSE, grid = FALSE, 
                                   gd_3 = graphData_conc,
                                   include_counts = FALSE)
  
  ggplot2::ggsave(toxPlot_wq, filename = file_to_save, width = 11, height = 9)

}


graphData_wq_tox <- function(graphData_tox, graphData_wq, file_to_save){
  

  graphData_tox$guide_side <- "ToxCast\nMaximum EAR per Site"
  graphData_wq$guide_side <- "Traditional\nMaximum Toxicity Quotient per Site"
  
  toxPlot_wq_tox <- combo_plot_matches(graphData_tox, graphData_wq, 
                                     thres_1 = NA, thres_2 = NA, 
                                     drop = TRUE, grid = FALSE,include_counts = FALSE)

  ggplot2::ggsave(toxPlot_wq_tox, filename = file_to_save, width = 11, height = 9)
  
}