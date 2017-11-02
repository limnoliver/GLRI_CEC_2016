library(toxEval)
library(dplyr)
library(ggplot2)

get_colors <- function(){
  cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e")
  return(cbValues)
}

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

bench_tox_conc_data <- function(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc, chem_info, benchmarks){
  
  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary$type <- "ToxCast"
  chemicalSummary_conc$type <- "Concentration"
  
  total_summary <- suppressWarnings(bind_rows(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc))
  
  chnm_df <- data.frame(CAS = chem_info$CAS, stringsAsFactors = FALSE) %>%
    left_join(distinct(select(ACC, CAS=casn, chnm)), by="CAS") %>%
    left_join(distinct(select(benchmarks, CAS, Compound)), by="CAS") %>%
    left_join(distinct(select(chem_info, CAS, `Chemical Name`)), by="CAS") %>%
    distinct()
  
  chnm_df$chnm[is.na(chnm_df$chnm)] <- chnm_df$Compound[is.na(chnm_df$chnm)]
  chnm_df$chnm[is.na(chnm_df$chnm)] <- chnm_df$`Chemical Name`[is.na(chnm_df$chnm)]
  
  graphData <-  total_summary %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark","Concentration"))) %>%
    left_join(chnm_df, by="CAS")
  
  orderClass_df <- toxEval:::orderClass(graphData)
  
  orderChem_tc <- toxEval:::orderChem(filter(graphData, type == "ToxCast"), orderClass_df)
  orderChem_bm <- toxEval:::orderChem(filter(graphData, type == "Benchmark"), orderClass_df) 
  orderChem_conc <- toxEval:::orderChem(filter(graphData, type == "Concentration"), orderClass_df) 

  orderedChems <- full_join(orderChem_tc, orderChem_bm, by=c("chnm","Class")) %>%
    full_join(orderChem_conc, by=c("chnm","Class")) %>%
    mutate(chnm = as.character(chnm)) %>%
    arrange(Class, median.x, median.y, median)
  
  graphData$chnm <- factor(graphData$chnm,
                           levels = orderedChems$chnm)
  graphData$Class <- factor(graphData$Class,
                            levels = orderClass_df$Class)
  return(graphData)
}

benchmark_tox_figure <- function(graphData, chemicalSummary, chemicalSummary_bench, cbValues, file_out){

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
    scale_y_log10(labels=fancyNumbers, breaks = c(1 %o% 10^(-8:1)))  +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=1) +
    facet_grid(. ~ type, scales = "free") +
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
  
  layout_stuff <- ggplot_build(toxPlot_All)
  layout_stuff$layout$panel_ranges[[1]]$x.range[1]
  
  astrictData_tox$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data = astrictData_bench, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)+
    geom_text(data = astrictData_tox, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)
  
  layout_stuff <- ggplot_build(toxPlot_All_withLabels)
  layout_stuff$layout$panel_ranges[[1]]$x.range[1]
  
  countNonZero$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  
  toxPlot_All_withLabels <- toxPlot_All_withLabels +
    geom_text(data=countNonZero, 
              aes(x= chnm, label = nonZero, y=y), size=1.75) 
  
  ggsave(toxPlot_All_withLabels, filename = file_out, width = 11, height = 9)
  
}

benchmark_tox_conc_figure <- function(graphData_all_3, chemicalSummary, chemicalSummary_bench, cbValues, file_out){

  countNonZero <- graphData_all_3 %>%
    select(CAS, chnm, Class, maxEAR, type) %>%
    group_by(CAS, chnm, Class, type) %>%
    summarize(nonZero = as.character(sum(maxEAR>0))) %>%
    ungroup() %>%
    select(-type) %>%
    distinct() %>%
    mutate(type = factor("ToxCast", levels = c("ToxCast","Benchmark","Concentration")),
           y=10^-8)
  
  astrictData_tox <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor("ToxCast", levels = c("ToxCast","Benchmark","Concentration"))) %>%
    filter(!(CAS %in% unique(chemicalSummary$CAS)))
  
  astrictData_bench <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor("Benchmark", levels = c("ToxCast","Benchmark"))) %>%
    filter(!(CAS %in% unique(chemicalSummary_bench$CAS)))
  
  toxPlot_All <- ggplot(data=graphData_all_3) +
    scale_y_log10(labels=fancyNumbers, breaks = c(1 %o% 10^(-8:1)))  +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=1) +
    facet_grid(. ~ type, scales = "free") +
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
          legend.key.height = unit(1,"line"),
          legend.position="bottom",
          legend.justification = "left") +
    scale_fill_manual(values = cbValues, drop=FALSE) 
  
  layout_stuff <- ggplot_build(toxPlot_All)
  
  astrictData_tox$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  astrictData_bench$y <- 10^(layout_stuff$layout$panel_ranges[[2]]$x.range[1])
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data = astrictData_bench, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)+
    geom_text(data = astrictData_tox, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)
  
  layout_stuff <- ggplot_build(toxPlot_All_withLabels)
  layout_stuff$layout$panel_ranges[[1]]$x.range[1]
  
  countNonZero$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  
  toxPlot_All_withLabels <- toxPlot_All_withLabels +
    geom_text(data=countNonZero, 
            aes(x= chnm, label = nonZero, y=y), size=1.75) 
  
  ggsave(toxPlot_All_withLabels, filename = file_out, width = 11, height = 9)
  
}

class_figures <- function(graphData_all_3, chemicalSummary, chemicalSummary_bench, cbValues, file_out){

  plotted <- data.frame(class = levels(graphData_all_3$Class),
                        file = "",stringsAsFactors = FALSE)
  
  for(i in seq_along(levels(graphData_all_3$Class))){
  
    plot_class <- levels(graphData_all_3$Class)[i]
    sub_data <- filter(graphData_all_3, Class == plot_class)
    
    countNonZero <- sub_data %>%
      select(CAS, chnm, Class, maxEAR, type) %>%
      group_by(CAS, chnm, Class, type) %>%
      summarize(nonZero = as.character(sum(maxEAR>0))) %>%
      ungroup() %>%
      select(-type) %>%
      distinct() %>%
      mutate(type = factor("ToxCast", levels = c("ToxCast","Benchmark","Concentration")),
             y=10^-8)
    
    astrictData_tox <- countNonZero %>%
      mutate(y = 10^-7.5,
             askt = "*",
             type = factor("ToxCast", levels = c("ToxCast","Benchmark","Concentration"))) %>%
      filter(!(CAS %in% unique(chemicalSummary$CAS)))
    
    astrictData_bench <- countNonZero %>%
      mutate(y = 10^-7.5,
             askt = "*",
             type = factor("Benchmark", levels = c("ToxCast","Benchmark"))) %>%
      filter(!(CAS %in% unique(chemicalSummary_bench$CAS)))
    
    toxPlot_All <- ggplot(data=sub_data) +
      scale_y_log10(labels=fancyNumbers, breaks = c(1 %o% 10^(-8:2)))  +
      geom_boxplot(aes(x=chnm, y=maxEAR),
                   lwd=0.1,outlier.size=1, fill = cbValues[i]) +
      facet_grid(. ~ type, scales = "free") +
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
      theme(legend.position="none") +
      ggtitle(label = plot_class)
    
    layout_stuff <- ggplot_build(toxPlot_All)
    
    if(nrow(astrictData_tox) > 0){
      astrictData_tox$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
    }
    
    if(nrow(astrictData_bench) > 0){
      astrictData_bench$y <- 10^(layout_stuff$layout$panel_ranges[[2]]$x.range[1])
    }
    
    
    toxPlot_All_withLabels <- toxPlot_All +
      geom_text(data = astrictData_bench, aes(x=chnm, label=askt, y=y),
                size=5, vjust = 0.70)+
      geom_text(data = astrictData_tox, aes(x=chnm, label=askt, y=y),
                size=5, vjust = 0.70)
    
    layout_stuff <- ggplot_build(toxPlot_All_withLabels)
    layout_stuff$layout$panel_ranges[[1]]$x.range[1]
    
    countNonZero$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
    
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=countNonZero, 
                aes(x= chnm, label = nonZero, y=y)) 
    
    file_to_save <- paste0("figure/Compare_",plot_class,".png")
    ggsave(toxPlot_All_withLabels, filename = file_to_save, width = 11, height = 9)
    plotted$file[i] <- file_to_save
  }
  
  write.csv(plotted, file = file_out, row.names = FALSE)
}