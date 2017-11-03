library(toxEval)
library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)

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

  bench_stuff <- distinct(select(benchmarks, CAS, Compound))
  bench_stuff <- bench_stuff[!duplicated(bench_stuff$CAS),]
  
  chnm_df <- data.frame(CAS = chem_info$CAS, stringsAsFactors = FALSE) %>%
    left_join(distinct(select(ACC, CAS=casn, chnm)), by="CAS") %>%
    left_join(bench_stuff, by="CAS") %>%
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
  
  bench_stuff <- distinct(select(benchmarks, CAS, Compound))
  bench_stuff <- bench_stuff[!duplicated(bench_stuff$CAS),]
  

  chnm_df <- data.frame(CAS = chem_info$CAS, stringsAsFactors = FALSE) %>%
    left_join(distinct(select(ACC, CAS=casn, chnm)), by="CAS") %>%
    left_join(bench_stuff, by="CAS") %>%
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

graph_tox_bench <- function(graphData, facet_levels, tox_chems, bench_chems, fill_boxes = "Dynamic"){
  
  countNonZero <- graphData %>%
    select(CAS, chnm, Class, maxEAR,type) %>%
    group_by(CAS, chnm, Class,type) %>%
    summarize(nonZero = sum(maxEAR>0)) %>%
    ungroup() %>%
    select(-type) %>%
    distinct() %>%
    mutate(type = factor(facet_levels[1], levels = facet_levels))

  astrictData_tox <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor(facet_levels[1], levels = facet_levels)) %>%
    filter(!(CAS %in% tox_chems))
  
  astrictData_bench <- countNonZero %>%
    mutate(y = 10^-7.5,
           askt = "*",
           type = factor(facet_levels[2], levels = facet_levels)) %>%
    filter(!(CAS %in% bench_chems))
  
  toxPlot_All <- ggplot(data=graphData) +
    scale_y_log10(labels=fancyNumbers, breaks = c(1 %o% 10^(-8:1))) 
  
  if(fill_boxes == "Dynamic"){
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=1) 
  } else {
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=maxEAR),
                   lwd=0.1,outlier.size=1, fill=fill_boxes) 
  }
  
  toxPlot_All <- toxPlot_All +
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
          legend.key.height = unit(1,"line")) 
  
  layout_stuff <- ggplot_build(toxPlot_All)
  layout_stuff$layout$panel_ranges[[1]]$x.range[1]
  
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
              aes(x= chnm, label = nonZero, y=y), size=1.75) 
  
  return(toxPlot_All_withLabels)
}

benchmark_tox_figure <- function(graphData, chemicalSummary, chemicalSummary_bench, cbValues, file_out){

  toxPlot_All_withLabels <- graph_tox_bench(graphData, 
                                            facet_levels = c("ToxCast","Benchmark"),
                                            tox_chems = unique(chemicalSummary$CAS),
                                            bench_chems = unique(chemicalSummary_bench$CAS),
                                            fill_boxes = "Dynamic") +
    scale_fill_manual(values = cbValues, drop=FALSE) 
  
  ggsave(toxPlot_All_withLabels, filename = file_out, width = 11, height = 9)
  
}

benchmark_tox_conc_figure <- function(graphData_all_3, chemicalSummary, chemicalSummary_bench, cbValues, file_out){

  toxPlot_All_withLabels <- graph_tox_bench(graphData_all_3, 
                                            facet_levels = c("ToxCast","Benchmark","Concentration"),
                                            tox_chems = unique(chemicalSummary$CAS),
                                            bench_chems = unique(chemicalSummary_bench$CAS),
                                            fill_boxes = "Dynamic") +
    scale_fill_manual(values = cbValues, drop=FALSE) +
    theme(legend.position="bottom",
          legend.justification = "left")
  
  ggsave(toxPlot_All_withLabels, filename = file_out, width = 11, height = 9)
  
}

class_figures <- function(graphData_all_3, chemicalSummary, chemicalSummary_bench, cbValues, file_out){

  plotted <- data.frame(class = levels(graphData_all_3$Class),
                        file = "",stringsAsFactors = FALSE)
  
  for(i in seq_along(levels(graphData_all_3$Class))){
  
    plot_class <- levels(graphData_all_3$Class)[i]
    sub_data <- filter(graphData_all_3, Class == plot_class)
    
    toxPlot_All_withLabels <- graph_tox_bench(sub_data, 
                                              facet_levels = c("ToxCast","Benchmark","Concentratino"),
                                              tox_chems = unique(chemicalSummary$CAS),
                                              bench_chems = unique(chemicalSummary_bench$CAS),
                                              fill_boxes = cbValues[i])
    
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      theme(legend.position="none") +
      ggtitle(label = plot_class)
    
    file_to_save <- paste0("figure/Compare_",plot_class,".png")
    ggsave(toxPlot_All_withLabels, filename = file_to_save, width = 11, height = 9)
    plotted$file[i] <- file_to_save
  }
  
  write.csv(plotted, file = file_out, row.names = FALSE)
}

plot_class_summaries <- function(file_out, chemicalSummary, category, title_words){

  if(category == "Chemical Class"){
    chemicalSummary$category <- chemicalSummary$Class
  } else {
    chemicalSummary$category <- chemicalSummary$Bio_category
  }
  
  tox_plot <- plot_tox_boxplots(chemicalSummary, category, mean_logic = FALSE)
  
  tox_plot <-tox_plot +
    ggtitle(label = title_words)
  
  ggsave(tox_plot, filename = file_out, width = 7, height = 7)
  
}

plot_facet_class <- function(target_name, chemicalSummary, chemicalSummary_bench, chemicalSummary_conc){

  chemicalSummary$type <- "ToxCast"
  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary_conc$type <- "Concentration"
  
  tots <- bind_rows(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc)
  
  tots$type <- factor(tots$type, levels = c("ToxCast","Benchmark","Concentration"))
  tots$Class <- factor(tots$Class, levels = rev(levels(chemicalSummary$Class)))
  
  tots$category <- tots$Class
  
  class_plot <- plot_tox_boxplots_facet(tots)
  
  class_plot <- class_plot +
    facet_grid(. ~ type, scales = "free") 
  
  ggsave(filename = target_name, plot = class_plot, width = 11, height = 5)
  
}

plot_tox_boxplots_facet <- function(chemicalSummary){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"

  graphData_df <- chemicalSummary %>%
    group_by(site,date, category, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, category, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() 
  
  bioPlot <- ggplot(data = graphData_df)+
    coord_flip() +
    theme_bw() +
    xlab("") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          axis.text.y = element_text(color = "black", vjust = 0.2), 
          axis.text.x = element_text(color = "black", vjust = 0, margin = margin(-0.5,0,0,0)))

  bioPlot <- bioPlot + 
    geom_boxplot(aes(x=category, y=maxEAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
    scale_y_log10("Maximum EAR Per Site",labels=fancyNumbers) 
  
  return(bioPlot)

}

plot_genes <- function(file_out, chem_info, chem_data, site_info, exclusions, AOP){
  
  ACClong <- get_ACC(chem_info$CAS)
  ACClong <- remove_flags(ACClong)
  
  cleaned_ep <- clean_endPoint_info(endPointInfo)
  filtered_ep <- filter_groups(cleaned_ep, groupCol = "intended_target_gene_symbol")

  chemicalSummary <- get_chemical_summary(ACClong,
                                          filtered_ep,
                                          chem_data, 
                                          site_info, 
                                          chem_info,
                                          exclusions)
  
  chemicalSummary <- select(chemicalSummary, -Class, -chnm) %>%
    left_join(AOP, by=c("Bio_category"="gene_symbol")) %>%
    rename(chnm = Bio_category,
           Class = AOP) %>%
    filter(!is.na(Class))
  
  x <- sapply(strsplit(chemicalSummary$Class, " leading to"), function(x)x[[1]])
  x <- sapply(strsplit(x, " leading to"), function(x)x[[1]])
  
  x <- gsub(" Beta-Oxidation Inhibition Leading to Steatosis, PPAR? activation","", x)
  x <- gsub(" induced by competitive antagonists of ionotropic GABA receptors","", x)
  x <- gsub("Peroxisomal Fatty Acid Beta-Oxidation Inhibition Leading to Steatosis, ","", x)
  x <- gsub(" leading neurotoxicity and excess acute toxicity","", x)
  x <- gsub("Upregulation of Thyroid Hormone Catabolism via Activation of Hepatic Nuclear Receptors, and Subsequent Adverse Neurodevelopmental Outcomes in Mammals , ","", x)
  x <- gsub(" function Subsequent to Estradiol Activation in the Fetal Testis","", x)
  x <- gsub(" function Subsequent to Estradiol Activation in the Fetal Testis","", x)
  
  chemicalSummary$Class <- x
  
  gd <- graph_chem_data(chemicalSummary)
  
  orderClass_df <- toxEval:::orderClass(gd)
  
  orderChem_df <- toxEval:::orderChem(gd, orderClass_df)

  chemicalSummary$chnm <- factor(chemicalSummary$chnm,
                                 levels = orderChem_df$chnm)    
  
  chemicalSummary$Class <- factor(chemicalSummary$Class,
                                  levels = rev(levels(orderChem_df$Class)))
  
  
  bioPlot <- plot_chemical_boxplots(chemicalSummary)
  
  ggsave(filename = file_out, plot = bioPlot, width = 11, height = 8)
}

plot_landuse <- function(target_name, sites){
  
  sites$shortName[is.na(sites$shortName)] <- "Saginaw"
  
  long_site <- select(sites, shortName, Urban, Ag=Ag..total, Forest, 
                      Wetland=Water..wetland, Other = Other.land.use) %>%
    gather(landuse, values, -shortName) 
  
  order_sites <- long_site %>%
    filter(landuse == "Ag") %>%
    arrange(desc(values))
  
  long_site$landuse <- factor(long_site$landuse, levels = c("Other","Wetland","Forest","Urban","Ag"))
  long_site$shortName <- factor(long_site$shortName, levels = order_sites$shortName)
  
  landuse <- ggplot(data = long_site, aes(x = shortName)) +
    geom_col(aes(y = values, fill = landuse)) +
    theme_minimal() +
    xlab("") +
    ylab("Land Use %") +
    theme(axis.text.x = element_blank())
  
  ggsave(filename = target_name, plot = landuse, width = 10, height = 7)
  
  return(landuse)
}

plot_stacks <- function(target_name, sites, chemicalSummary, site_info){
  
  sites$shortName[is.na(sites$shortName)] <- "Saginaw"
  
  long_site <- select(sites, shortName, Urban, Ag=Ag..total, Forest, 
                      Wetland=Water..wetland, Other = Other.land.use) %>%
    gather(landuse, values, -shortName) 
  
  order_sites <- long_site %>%
    filter(landuse == "Ag") %>%
    arrange(desc(values))
  
  long_site$landuse <- factor(long_site$landuse, levels = c("Other","Wetland","Forest","Urban","Ag"))
  long_site$shortName <- factor(long_site$shortName, levels = order_sites$shortName)
  
  chem_site$`Short Name` <- factor(chem_site$`Short Name`, levels = order_sites$shortName)
  
  stack_fig <- plot_tox_stacks(chemicalSummary, chem_site)
  
  ggsave(filename = target_name, plot = stack_fig, width = 10, height = 7)
  
  return(stack_fig)
}

combine_figures <- function(target_name, landuse_fig, stacK_fig){

  png(file = target_name, height = 900, width = 1100)
  grid.newpage()
  grid.draw(rbind(ggplotGrob(stacK_fig), ggplotGrob(landuse_fig),  size = "last")) 
  dev.off()

}
