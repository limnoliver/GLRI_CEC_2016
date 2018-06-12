endpoint_graph <- function(file_path, chemicalSummary, AOP, threshold, siteThreshold){

  library(dplyr)
  
  endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
    group_by(endPoint,site,date) %>%
    summarize(EARsum = sum(EAR)) %>%
    group_by(site,endPoint) %>%
    summarize(EARmax = max(EARsum)) %>%
    filter(EARmax >= threshold) %>%
    group_by(endPoint) %>%
    summarize(numSites = n_distinct(site)) %>%
    arrange(desc(numSites)) %>%
    filter(numSites >= siteThreshold)
  
  priority_endpoints <- endpoints_sites_hits$endPoint
  
  chemicalSummaryPriority <- filter(chemicalSummary, endPoint %in% priority_endpoints)
  
  eps_with_ids <- unique(AOP$endPoint)
  
  chemicalSummaryPriority$has_AOP <- "AOP Undefined"
  chemicalSummaryPriority$has_AOP[chemicalSummaryPriority$endPoint %in% eps_with_ids] <- "AOP Associated"
  
  endpointPlot <- plot_tox_endpoints_manuscript(chemicalSummaryPriority)
  
  gb <- ggplot2::ggplot_build(endpointPlot)
  gt <- ggplot2::ggplot_gtable(gb)
  
  gt$layout$clip[gt$layout$name=="panel"] <- "off"
  
  dir.create(file.path("plots"), showWarnings = FALSE)
  png(file_path, width = 1000, height = 800, res = 142)
  grid::grid.draw(gt)
  dev.off()
  
}