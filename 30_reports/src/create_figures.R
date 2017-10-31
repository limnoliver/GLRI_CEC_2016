library(toxEval)
library(dplyr)
library(ggplot2)


benchmark_tox <- function(chemicalSummary, 
                          chem_data, site_info, chem_info, 
                          benchmarks, exclusions, file_out){
  
  benchmarks <- benchmarks %>%
    rename(chnm = Compound,
           ACC_value = value) %>%
    filter(!is.na(CAS))
  
  filtered_ep <- select(benchmarks, endPoint) %>%
    distinct() %>%
    mutate(groupCol = "Aquatic Benchmark")
  
  chemicalSummary_bench <- get_chemical_summary(benchmarks,
                                          filtered_ep,
                                          chem_data, 
                                          site_info, 
                                          chem_info,
                                          exclusions)
  
  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary$type <- "ToxCast"
  
  total_summary <- suppressWarnings(bind_rows(chemicalSummary, chemicalSummary_bench))

  total_summary$Class <- factor(total_summary$Class, 
                                levels = levels(chemicalSummary$Class))  
  
  chnm_df <- data.frame(CAS = chem_info$CAS, stringsAsFactors = FALSE) %>%
    left_join(distinct(select(ACC, CAS=casn, chnm)), by="CAS")
  
  total_GD <- total_summary %>%
    select(-Bio_category, -shortName, -chnm) %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    rename(SiteID = site) %>%
    left_join(chnm_df, by="CAS") %>%
    left_join(distinct(select(benchmarks, CAS, ChemName=chnm)), by="CAS")
  
  total_GD$chnm[is.na(total_GD$chnm)] <- total_GD$ChemName[is.na(total_GD$chnm)]
  
  total_GD$chnm <- factor(total_GD$chnm, 
                          levels = c(levels(chemicalSummary$chnm),
                                     unique(total_GD$chnm)[!unique(total_GD$chnm) %in% levels(chemicalSummary$chnm)]))
  
  cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e")
  
  toxPlot_All <- ggplot(data=total_GD) +
    scale_y_log10(labels=fancyNumbers)  +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=1) +
    facet_grid(. ~ type, scales = "free", space = "free") +
    theme_bw() +
    scale_x_discrete(drop=TRUE) +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=7),
          axis.title=element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank()) +
    guides(fill=guide_legend(ncol=6)) +
    theme(legend.position="bottom",
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line")) +
    scale_fill_manual(values = cbValues, drop=FALSE) 
  
  ggsave(toxPlot_All, filename = file_out, width = 11, height = 9)
  
}