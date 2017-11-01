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
  
  graphData <-  total_summary %>%
    group_by(site,date,chnm, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, chnm, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark")))
  
  orderClass_df <- toxEval:::orderClass(graphData)
  
  orderChem_df <- toxEval:::orderChem(graphData, orderClass_df)
  
  graphData$chnm <- factor(graphData$chnm,
                           levels = orderChem_df$chnm)    
  
  cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e")
  
  countNonZero <- graphData %>%
    select(chnm, Class, maxEAR, type) %>%
    group_by(chnm, Class, type) %>%
    summarize(nonZero = as.character(sum(maxEAR>0))) %>%
    ungroup() %>%
    select(-type) %>%
    distinct() %>%
    mutate(type = factor("ToxCast", levels = c("ToxCast","Benchmark")),
           y=10^-8)

  astrictData_tox <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor("ToxCast", levels = c("ToxCast","Benchmark"))) %>%
    filter(!(chnm %in% unique(chemicalSummary$chnm)))
  
  astrictData_bench <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor("Benchmark", levels = c("ToxCast","Benchmark"))) %>%
    filter(!(chnm %in% unique(chemicalSummary_bench$chnm)))
    
  
  toxPlot_All <- ggplot(data=graphData) +
    scale_y_log10(labels=fancyNumbers, breaks = c(1 %o% 10^(-8:0)))  +
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
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, 
              aes(x= chnm, label = nonZero, y=y), size=1.75) +
    geom_text(data = astrictData_bench, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)+
    geom_text(data = astrictData_tox, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)
  
  ggsave(toxPlot_All_withLabels, filename = file_out, width = 11, height = 9)
  
}