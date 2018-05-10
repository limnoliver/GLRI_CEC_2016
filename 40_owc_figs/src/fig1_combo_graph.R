graphData_combos <- function(cs_1, cs_2, cs_3, file_to_save){
  
  graphData_tox <- graph_chem_data(cs_1)
  graphData_tox$guide_side <- "ToxCast\nMaximum EAR per Site"
  
  graphData_wq <- graph_chem_data(cs_2, sum_logic = FALSE)
  graphData_wq$guide_side <- "Traditional\nMaximum Toxicity Quotient per Site"
  
  graphData_conc <- graph_chem_data(cs_3, sum_logic = FALSE)
  graphData_conc$guide_side <- "Concentration [ug/L]"

  browser()
  
  toxPlot_wq <- combo_plot_matches(graphData_tox, graphData_wq, 
                                   thres_1 = NA, thres_2 = NA, 
                                   drop = FALSE, grid = FALSE, 
                                   gd_3 = graphData_conc,
                                   include_counts = FALSE)
  
  toxPlot_wq_tox <- combo_plot_matches(graphData_tox, graphData_wq, 
                                   thres_1 = NA, thres_2 = NA, 
                                   drop = TRUE, grid = FALSE)
  
  
  
}
