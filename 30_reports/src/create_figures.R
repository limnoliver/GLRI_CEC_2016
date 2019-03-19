library(toxEval)
library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
library(RColorBrewer)

get_colors <- function(chemicalSummary_conc){
  
  cbValues <- brewer.pal(n = length(levels(chemicalSummary_conc$Class)), name = 'Set1')
  # cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
  #               "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
  #               "#FFA500","#F4426e")

  names(cbValues) <- levels(chemicalSummary_conc$Class)
  return(cbValues)
}

bench_tox_data <- function(chemicalSummary, chemicalSummary_bench, chem_info, benchmarks){

  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary$type <- "ToxCast"
  
  total_summary <- suppressWarnings(bind_rows(chemicalSummary, chemicalSummary_bench))

  bench_stuff <- distinct(select(benchmarks, CAS, Compound))
  bench_stuff <- bench_stuff[!duplicated(bench_stuff$CAS),]
  
  chnm_df <- data.frame(CAS = chem_info$CAS, chnm = chem_info$`Chemical Name`,
                        stringsAsFactors = FALSE) 
  
  graphData_df <-  total_summary %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark"))) %>%
    left_join(chnm_df, by="CAS")
  
  orderClass_df <- toxEval:::orderClass(filter(graphData_df, type == "ToxCast"))
  
  orderChem_tc <- toxEval:::orderChem(filter(graphData_df, type == "ToxCast"), orderClass_df)
  orderChem_bm <- toxEval:::orderChem(filter(graphData_df, type == "Benchmark"), orderClass_df) 
  
  orderedChems <- full_join(orderChem_tc, orderChem_bm, by=c("chnm","Class")) %>%
    mutate(chnm = as.character(chnm)) %>%
    arrange(Class, median.x, median.y)
  
  graphData_df$chnm <- factor(graphData_df$chnm,
                           levels = orderedChems$chnm)
  graphData_df$Class <- factor(graphData_df$Class,
                            levels = orderClass_df$Class)
  return(graphData_df)
}

bench_tox_conc_data <- function(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc, chem_info, benchmarks){
  
  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary$type <- "ToxCast"
  chemicalSummary_conc$type <- "Concentration"
  
  total_summary <- suppressWarnings(bind_rows(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc))
  
  bench_stuff <- distinct(select(benchmarks, CAS, Compound))
  bench_stuff <- bench_stuff[!duplicated(bench_stuff$CAS),]
  
  chnm_df <- data.frame(CAS = chem_info$CAS, chnm = chem_info$`Chemical Name`,
                        stringsAsFactors = FALSE) 

  graphData_df <-  total_summary %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark","Concentration"))) %>%
    left_join(chnm_df, by="CAS")
  
  orderClass_df <- toxEval:::orderClass(filter(graphData_df, type == "ToxCast"))
  
  orderChem_tc <- toxEval:::orderChem(filter(graphData_df, type == "ToxCast"), orderClass_df)
  orderChem_bm <- toxEval:::orderChem(filter(graphData_df, type == "Benchmark"), orderClass_df) 
  orderChem_conc <- toxEval:::orderChem(filter(graphData_df, type == "Concentration"), orderClass_df) 

  orderedChems <- full_join(orderChem_tc, orderChem_bm, by=c("chnm","Class")) %>%
    full_join(orderChem_conc, by=c("chnm","Class")) %>%
    mutate(chnm = as.character(chnm)) %>%
    arrange(Class, median.x, median.y, median)
  
  graphData_df$chnm <- factor(graphData_df$chnm,
                           levels = orderedChems$chnm)
  graphData_df$Class <- factor(graphData_df$Class,
                            levels = orderClass_df$Class)
  return(graphData_df)
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
    scale_y_log10(labels=toxEval:::fancyNumbers, breaks = c(1 %o% 10^(-8:1))) 
  
  if(fill_boxes == "Dynamic"){
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 outlier.size=0.75)#, position = position_dodge2(preserve = "total")) 
  } else {
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=maxEAR),
                   lwd=0.1,outlier.size=0.5, fill=fill_boxes)#, position = position_dodge2(preserve = "total")) 
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
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank()) +
    theme(legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line")) 
  
  layout_stuff <- ggplot_build(toxPlot_All) 

  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    y_tox <- 10^(layout_stuff$layout$panel_scales_y[[1]]$range$range[1])
    y_bench <- 10^(layout_stuff$layout$panel_scales_y[[2]]$range$range[1])
  } else {
    y_tox <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
    y_bench <- 10^(layout_stuff$layout$panel_ranges[[2]]$x.range[1])
  } 

  if(nrow(astrictData_tox) > 0){
    astrictData_tox$y <- y_tox
  }
  if(nrow(astrictData_bench) > 0){
    astrictData_bench$y <- y_bench
  }
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data = astrictData_bench, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)+
    geom_text(data = astrictData_tox, aes(x=chnm, label=askt, y=y),
              size=5, vjust = 0.70)
  
  layout_stuff <- ggplot_build(toxPlot_All_withLabels)
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    nonZ_y <- 10^(layout_stuff$layout$panel_scales_y[[1]]$range$range[1])
  } else {
    nonZ_y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  } 
  
  countNonZero$y <- nonZ_y
  
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
    geom_boxplot(aes(x=category, y=maxEAR),lwd=0.1,outlier.size=0.75, fill = "steelblue") +
    scale_y_log10("Maximum EAR Per Site",labels=toxEval:::fancyNumbers) 
  
  return(bioPlot)

}

plot_genes <- function(file_out, chem_info, chem_data, site_info, exclusions, AOP){
  
  ACClong <- get_ACC(chem_info$CAS)
  ACClong <- remove_flags(ACClong)
  
  cleaned_ep <- clean_endPoint_info(end_point_info)
  filtered_ep <- filter_groups(cleaned_ep, groupCol = "intended_target_gene_symbol")

  chemicalSummary <- get_chemical_summary(tox_list = NULL,
                                          ACClong = ACClong,
                                          filtered_ep = filtered_ep,
                                          chem.data = chem_data, 
                                          chem.site = site_info,
                                          chem.info = chem_info,
                                          exclusion = exclusions)
  
  chemicalSummary <- select(chemicalSummary, -Class, -chnm) %>%
    left_join(AOP, by=c("Bio_category"="gene_symbol")) %>%
    rename(chnm = Bio_category,
           Class = AOP) 
  
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
  chemicalSummary$Class[is.na(chemicalSummary$Class)] <- "Not Classified"

  chemicalSummary$chnm <- factor(chemicalSummary$chnm,
                                 levels = orderChem_df$chnm)    
  
  chemicalSummary$Class <- factor(chemicalSummary$Class,
                                  levels = c(rev(levels(orderChem_df$Class)),"Not Classified"))
  
  bioPlot <- plot_chemical_boxplots(chemicalSummary) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    guides(fill=guide_legend(ncol=4)) 
  
  ggsave(filename = file_out, plot = bioPlot, width = 11, height = 20)
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
    theme(axis.text.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  
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
  
  site_info$`Short Name` <- factor(site_info$`Short Name`, levels = order_sites$shortName)
  
  stack_fig <- plot_tox_stacks(chemicalSummary, site_info)
  
  ggsave(filename = target_name, plot = stack_fig, width = 10, height = 7)
  
  return(stack_fig)
}

combine_figures <- function(target_name, landuse_fig, stacK_fig){

  png(file = target_name, height = 900, width = 1100)
  grid.newpage()
  grid.draw(rbind(ggplotGrob(stacK_fig), ggplotGrob(landuse_fig),  size = "last")) 
  dev.off()

}

conc_plot <- function(chemicalSummary_conc, target_name){
  
  chem_bp <- plot_chemical_boxplots(chemicalSummary_conc, plot_ND = FALSE)
  
  chem_bp <- chem_bp +
    theme(legend.position = "right",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(color = "black"),
          legend.text=element_text(size=14),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    guides(fill=guide_legend(ncol=1)) 
  
  chem_bp$layers <- list(chem_bp$layers[[1]])
  
  ggsave(filename = target_name, plot = chem_bp, width = 7, height = 5)
  
}

plot_two <- function(graphData_b_c, chemicalSummary_bench, chemicalSummary_conc, cbValues, facet_labels){
  
  graphData_b_c <- filter(graphData_b_c, CAS %in% unique(chemicalSummary_bench$CAS))
  
  graphData_b_c$type <- factor(graphData_b_c$type, labels = facet_labels)
  
  countNonZero <- graphData_b_c %>%
    select(CAS, chnm, Class, maxEAR,type) %>%
    group_by(CAS, chnm, Class,type) %>%
    summarize(nonZero = sum(maxEAR>0)) %>%
    ungroup() %>%
    filter(type == levels(graphData_b_c$type)[1])
  
  toxPlot_All <- ggplot(data=graphData_b_c) +
    scale_y_log10(labels=toxEval:::fancyNumbers, breaks = c(1 %o% 10^(-8:1))) +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=0.75)+ #, position = position_dodge2(preserve = "total"))  +
    facet_grid(. ~ type, scales = "free") +
    theme_bw() +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = cbValues, drop=TRUE) +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=7),
          axis.title=element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank()) +
    theme(legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line")) 
  
  layout_stuff <- ggplot_build(toxPlot_All)
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    countNonZero$y_tox <- 10^(layout_stuff$layout$panel_scales_y[[1]]$range$range[1])
    
  } else {
    countNonZero$y_tox <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  } 
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, 
              aes(x= chnm, label = nonZero, y=y_tox), size=1.75) 
  
  return(toxPlot_All_withLabels)
}

benchmark_conc_figure <- function(graphData_b_c, chemicalSummary_bench, chemicalSummary_conc, cbValues, file_out){
  
  facet_labels <- c("Benchmark","Concentration")

  toxPlot_All <- plot_two(graphData_b_c, chemicalSummary_bench, chemicalSummary_conc, cbValues, facet_labels)
  
  ggsave(toxPlot_All, filename = file_out, width = 10, height = 5)
  
}

tox_conc_figure <- function(graphData_b_c, chemicalSummary_bench, chemicalSummary_conc, cbValues, file_out){
  
  facet_labels <- c("ToxCast","Concentration")
  toxPlot_All <- plot_two(graphData_b_c, chemicalSummary_bench, chemicalSummary_conc, cbValues, facet_labels)
  
  ggsave(toxPlot_All, filename = file_out, width = 10, height = 5)
  
}

all_3_cleaned <- function(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc, chem_info, cbValues, target_name){

  chemicalSummary$type <- "ToxCast"
  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary_conc$type <- "Concentration"
  
  chems_to_keep <- unique(c(chemicalSummary$CAS,chemicalSummary_bench$CAS))
  
  chemicalSummary_bench <- filter(chemicalSummary_bench, CAS %in% chems_to_keep)
  chemicalSummary_conc <- filter(chemicalSummary_conc, CAS %in% chems_to_keep)

  all_data <- rbind(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc)

  chnm_df <- select(chem_info, CAS, chnm=`Chemical Name`)
  # Order by tox:

  graphData_df <-  all_data %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark","Concentration"))) %>%
    left_join(chnm_df, by="CAS")

  orderClass_df <- toxEval:::orderClass(filter(graphData_df, type == "ToxCast"))

  orderChem_tc <- toxEval:::orderChem(filter(graphData_df, type == "ToxCast"), orderClass_df)
  orderChem_bm <- toxEval:::orderChem(filter(graphData_df, type == "Benchmark"), orderClass_df) 
  orderChem_conc <- toxEval:::orderChem(filter(graphData_df, type == "Concentration"), orderClass_df) 
  
  orderedChems <- full_join(orderChem_tc, orderChem_bm, by=c("chnm","Class")) %>%
    mutate(chnm = as.character(chnm)) %>%
    arrange(Class, median.x, median.y)

  graphData_df$chnm <- factor(graphData_df$chnm,
                           levels = orderedChems$chnm)
  graphData_df$Class <- factor(graphData_df$Class,
                            levels = orderClass_df$Class)
  
  astrictData_tox <- select(graphData_df, chnm, CAS) %>%
    distinct() %>%
    mutate(askt = "*",
           type = factor("ToxCast", levels = levels(graphData_df$type))) %>%
    filter(!(CAS %in% unique(chemicalSummary$CAS)))
  
  astrictData_bench <-select(graphData_df, chnm, CAS) %>%
    distinct() %>%
    mutate(askt = "*",
           type = factor("Benchmark", levels = levels(graphData_df$type))) %>%
    filter(!(CAS %in% unique(chemicalSummary_bench$CAS)))

  toxPlot_All <- ggplot(data=graphData_df) +
    scale_y_log10(labels=toxEval:::fancyNumbers, breaks = c(1 %o% 10^(-8:1))) +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class, color = Class),
                 outlier.size=0.75) + #, position = position_dodge2(preserve = "total")) +
    facet_grid(. ~ type, scales = "free") +
    theme_bw() +
    scale_x_discrete(drop=TRUE) +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=7),
          axis.title=element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line")) +
    scale_fill_manual(values = cbValues, drop=TRUE) +
    scale_color_manual(values = cbValues, drop = TRUE) +
    theme(legend.position = "right",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(color = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.text=element_text(size=14)) +
    guides(fill=guide_legend(ncol=1))
  
  layout_stuff <- ggplot_build(toxPlot_All)
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    astrictData_tox$y <- 10^(layout_stuff$layout$panel_scales_y[[1]]$range$range[1])
    astrictData_bench$y <- 10^(layout_stuff$layout$panel_scales_y[[2]]$range$range[1])
  } else {
    astrictData_tox$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
    astrictData_bench$y <- 10^(layout_stuff$layout$panel_ranges[[2]]$x.range[1])
  } 
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=astrictData_tox, 
              aes(x= chnm, label = askt, y=y)) +
    geom_text(data=astrictData_bench, 
              aes(x= chnm, label = askt, y=y))
  
  ggsave(toxPlot_All_withLabels, filename = target_name, width = 10, height = 5)
  
}

all_3_filtered <- function(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc, chem_info, cbValues, target_name){
  
  chemicalSummary$type <- "ToxCast"
  chemicalSummary_bench$type <- "Benchmark"
  chemicalSummary_conc$type <- "Concentration"
  
  chems_to_keep <- unique(c(chemicalSummary$CAS,chemicalSummary_bench$CAS))
  
  chemicalSummary_bench <- filter(chemicalSummary_bench, CAS %in% chems_to_keep)
  chemicalSummary_conc <- filter(chemicalSummary_conc, CAS %in% chems_to_keep)
  
  all_data <- rbind(chemicalSummary, chemicalSummary_bench)

  graphData_df <-  all_data %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame()
  
  median_stuff <- graphData_df %>%
    group_by(CAS, type) %>%
    summarise(med = median(maxEAR)) %>%
    filter(med > 0.001) %>%
    select(-type) %>%
    arrange(med)

  CAS_to_use <- median_stuff$CAS[!duplicated(median_stuff$CAS)]
#######################################################  
  all_data <- rbind(chemicalSummary, chemicalSummary_bench, chemicalSummary_conc)
  
  graphData_df <-  all_data %>%
    filter(CAS %in% CAS_to_use) %>%
    group_by(site,date,CAS, Class, type) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class, type) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor(type, levels = c("ToxCast","Benchmark","Concentration")))  %>%
    left_join(select(chem_info, CAS, chnm = `Chemical Name`), by="CAS")
  
  orderClass_df <- toxEval:::orderClass(filter(graphData_df, type == "ToxCast"))
  
  orderChem_tc <- toxEval:::orderChem(filter(graphData_df, type == "ToxCast"), orderClass_df)
  orderChem_bm <- toxEval:::orderChem(filter(graphData_df, type == "Benchmark"), orderClass_df) 
  orderChem_conc <- toxEval:::orderChem(filter(graphData_df, type == "Concentration"), orderClass_df) 
  
  orderedChems <- full_join(orderChem_tc, orderChem_bm, by=c("chnm","Class")) %>%
    mutate(chnm = as.character(chnm)) %>%
    arrange(Class, median.x, median.y)
  
  graphData_df$chnm <- factor(graphData_df$chnm,
                              levels = orderedChems$chnm)
  graphData_df$Class <- factor(graphData_df$Class,
                               levels = orderClass_df$Class)
  
  astrictData_tox <- select(graphData_df, chnm, CAS) %>%
    distinct() %>%
    mutate(askt = "*",
           type = factor("ToxCast", levels = levels(graphData_df$type))) %>%
    filter(!(CAS %in% unique(chemicalSummary$CAS)))

  dummy2 <- data.frame(type = factor(c("ToxCast", "Benchmark"),
                                     levels = levels(graphData_df$type)),
                       thresh = 10^-3)
  
  toxPlot_All <- ggplot(data=graphData_df) +
    scale_y_log10(labels=toxEval:::fancyNumbers, breaks = c(1 %o% 10^(-8:1))) +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class, color = Class),
                 outlier.size=0.75) +#, position = position_dodge2(preserve = "total")) +
    geom_hline(data = dummy2, aes(yintercept = thresh), linetype = "longdash") +
    facet_grid(. ~ type, scales = "free") +
    theme_bw() +
    scale_x_discrete(drop=TRUE) +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=7),
          axis.title=element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line")) +
    scale_fill_manual(values = cbValues, drop=TRUE) +
    scale_color_manual(values = cbValues, drop = TRUE) +
    theme(legend.position = "right",
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(color = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.text=element_text(size=14)) +
    guides(fill=guide_legend(ncol=1))

  astrictData_tox$y <- 10^-6
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=astrictData_tox, 
              aes(x= chnm, label = askt, y=y)) 
  
  layout_stuff <- ggplot_build(toxPlot_All_withLabels)
  
  countNonZero <- chemicalSummary_conc %>%
    filter(CAS %in% CAS_to_use) %>%
    group_by(site,date,CAS, Class) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, CAS, Class) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() %>%
    mutate(type = factor("ToxCast", levels = c("ToxCast","Benchmark","Concentration")))  %>%
    left_join(select(chem_info, CAS, chnm = `Chemical Name`), by="CAS") %>%
    select(CAS, chnm, Class, maxEAR,type) %>%
    group_by(CAS, chnm, Class,type) %>%
    summarize(nonZero = sum(maxEAR>0)) %>%
    ungroup() 
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    countNonZero$y <- 10^(layout_stuff$layout$panel_scales_y[[1]]$range$range[1])
  } else {
    countNonZero$y <- 10^(layout_stuff$layout$panel_ranges[[1]]$x.range[1])
  } 
  
  toxPlot_All_withLabels <- toxPlot_All_withLabels +
    geom_text(data=countNonZero, 
              aes(x= chnm, label = nonZero, y=y), size=2)
  
  ggsave(toxPlot_All_withLabels, filename = target_name, width = 10, height = 5)
  
  
}