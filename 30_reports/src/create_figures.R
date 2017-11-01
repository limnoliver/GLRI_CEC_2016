library(toxEval)
library(dplyr)
library(ggplot2)


bench_tox_data <- function(chemicalSummary, chemicalSummary_bench, chem_info, benchmarks){

  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary$type <- "ToxCast"
  
  total_summary <- suppressWarnings(bind_rows(chemicalSummary, chemicalSummary_bench))
  
  chnm_df <- data.frame(CAS = chem_info$CAS, stringsAsFactors = FALSE) %>%
    left_join(distinct(select(ACC, CAS=casn, chnm)), by="CAS") %>%
    left_join(distinct(select(benchmarks, CAS, Compound))) %>%
    distinct()
  
  chnm_df$chnm[is.na(chnm_df$chnm)] <- chnm_df$Compound[is.na(chnm_df$chnm)]
  
  graphData <-  total_summary %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark"))) %>%
    left_join(chnm_df, by="CAS")
  
  orderClass_df <- toxEval:::orderClass(graphData)
  
  orderChem_tc <- toxEval:::orderChem(filter(graphData, type == "ToxCast"), orderClass_df)
  orderChem_bm <- toxEval:::orderChem(filter(graphData, type == "Benchmark"), orderClass_df) 
  
  orderedChems <- full_join(orderChem_tc, orderChem_bm, by=c("chnm","Class")) %>%
    mutate(chnm = as.character(chnm)) %>%
    arrange(Class, median.x, median.y)
  
  graphData$chnm <- factor(graphData$chnm,
                           levels = orderedChems$chnm)
  graphData$Class <- factor(graphData$Class,
                            levels = orderClass_df$Class)
  return(graphData)
}

benchmark_tox <- function(graphData, chemicalSummary, chemicalSummary_bench, file_out){
  
  cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e")
  
  countNonZero <- graphData %>%
    select(CAS, chnm, Class, maxEAR, type) %>%
    group_by(CAS, chnm, Class, type) %>%
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
    filter(!(CAS %in% unique(chemicalSummary$CAS)))
  
  astrictData_bench <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor("Benchmark", levels = c("ToxCast","Benchmark"))) %>%
    filter(!(CAS %in% unique(chemicalSummary_bench$CAS)))
    
  
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
    theme(legend.background = element_rect(fill = "transparent", colour = "transparent"),
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