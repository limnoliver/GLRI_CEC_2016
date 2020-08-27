extract_gene_info <- function(file_in) {
  load(file_in)
  
  return(gene_info)
}

merge_with_aops <- function(top_ends, genes, top_mix, sites, site_threshold, sample_threshold, aop){
  
  filt_top_ends <- filter(top_ends, n_sites > site_threshold & n_hits > sample_threshold)  %>%
    left_join(select(aop, -endPointID, -AOPID, -KeyID))
  
  sites <- sites %>%
    mutate(disturbance = Ag..total + Urban) %>%
    select(site = site_no, 
           shortName, 
           dominant_lu = Dominant.land.use.,
           disturbance)
  
  
  browser()
  lu_top_ends <- top_mix %>%
    group_by(site, endPoint) %>%
    summarize(hits = n(),
              max_EARmix = max(sum_ear_endpoint)) %>%
    left_join(sites)
  
  
  
  # by_gene <- top_mix %>%
  #   filter(sum_ear_endpoint > 0.01) %>%
  #   left_join(genes) %>%
  #   group_by(site, geneSymbol) %>%
  #   summarize(n_months = length(unique(lubridate::month(date)))) %>%
  #   left_join(sites) %>%
  #   filter(!is.na(geneSymbol))

  
  tops <- filter(top_mix, sum_ear_endpoint > 0.01)
  test <- filter(genes, endPoint %in% unique(tops$endPoint))
  
  test_group <- test %>%
    group_by(geneSymbol) %>%
    summarize(n_endpoints = n())
    
  
  by_endpoint <- top_mix %>%
    group_by(endPoint) %>%
    mutate(n = n()) %>% ungroup() %>%
    filter(sum_ear_endpoint > 0.01) %>%
    group_by(site, endPoint) %>%
    summarize(n_months = length(unique(lubridate::month(date)))) %>%
    left_join(sites)
  
  site_order <- sites %>%
    filter(site %in% unique(by_endpoint$site)) %>%
    arrange(disturbance) %>%
    pull(shortName)
  
  #by_gene$shortName <- factor(by_gene$shortName, levels = site_order)
  #by_gene$dominant_lu <- factor(by_gene$dominant_lu, levels = c('Urban', 'Crops', 'AgMix'))
  
  by_endpoint$shortName <- factor(by_endpoint$shortName, levels = site_order)
  by_endpoint$dominant_lu <- factor(by_endpoint$dominant_lu, levels = c('Urban', 'Crops', 'AgMix'))
  
  # ggplot(by_gene, aes(y = shortName, x = geneSymbol)) +
  #   geom_tile(aes(fill = n_months), na.rm = FALSE) +
  #   facet_wrap(vars(dominant_lu), nrow = 4, scales = 'free_y', strip.position = 'right') +
  #   scale_fill_viridis_c(direction = -1, na.value = 'white') +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #   labs(x = 'Gene symbol', y = '', fill = 'Months with \nexceedances')
  
  p_end <- ggplot(by_endpoint, aes(y = shortName, x = endPoint)) +
    geom_tile(aes(fill = n_months), na.rm = FALSE) +
    facet_wrap(vars(dominant_lu), nrow = 4, scales = 'free_y', strip.position = 'right') +
    scale_fill_viridis_c(direction = -1, na.value = 'white') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = 'ToxCast assay', y = '', fill = 'Months with \nexceedances')

  ggsave('figures/ms_supplement_figures/top_endpoints_across_sites.png', p_end, height = 6, width = 9)

  # create top endpoints by land use
  by_endpoint <- group_by(by_endpoint, dominant_lu, endPoint) %>%
    summarize(n_hits = sum(n_months))
  
  # bio evaluation
  lu_top_ends <- filter(by_endpoint, n_hits >= 10)
  lu_top_ends_agmix <- filter(by_endpoint, n_hits >=4)
  # now find mixtures that were in these hits
  head(top_mix)
  chem_mixes <- top_mix %>%
    group_by(endPoint) %>%
    mutate(n = n()) %>% ungroup() %>%
    filter(sum_ear_endpoint > 0.01) %>%
    left_join(select(sites, site, dominant_lu), by = 'site')
  
  urban <- filter(chem_mixes) %>%
    filter(dominant_lu %in% 'Urban') %>%
    filter(endPoint %in% lu_top_ends$endPoint[lu_top_ends$dominant_lu == 'Urban']) %>%
    group_by(endPoint) %>%
    summarize(c_parents = paste(unique(unlist(strsplit(contr_parents, ', '))), collapse = ', '),
              c_all_chems = paste(unique(unlist(strsplit(contr_chems, ', '))), collapse = ', '),
              n = n())
  agmix <- filter(chem_mixes) %>%
    filter(dominant_lu %in% 'AgMix') %>%
    filter(endPoint %in% lu_top_ends_agmix$endPoint[lu_top_ends_agmix$dominant_lu == 'AgMix']) %>%
    group_by(endPoint) %>%
    summarize(c_parents = paste(unique(unlist(strsplit(contr_parents, ', '))), collapse = ', '),
              c_all_chems = paste(unique(unlist(strsplit(contr_chems, ', '))), collapse = ', '),
              n = n())
  crops <- filter(chem_mixes) %>%
    filter(dominant_lu %in% 'Crops') %>%
    filter(endPoint %in% lu_top_ends$endPoint[lu_top_ends$dominant_lu == 'Crops']) %>%
    group_by(endPoint) %>%
    summarize(c_parents = paste(unique(unlist(strsplit(contr_parents, ', '))), collapse = ', '),
              c_all_chems = paste(unique(unlist(strsplit(contr_chems, ', '))), collapse = ', '),
              n = n())
  
  
  lu_worst <- lu_top_ends %>%
    group_by(site, shortName, dominant_lu) %>%
    mutate(max_hits = max(hits)) %>% 
    ungroup() %>%
    filter(hits == max_hits)
  

  
  head(top_mix)
  
  test <- left_join(lu_worst, distinct(select(top_mix, site, endPoint, n_contr_chems, n_contr_parents, contr_chems, contr_parents)))
  
  lu_worst_mix <-  lu_top_ends %>%
    filter(n_contr_parents > 1) %>%
    group_by(site, shortName, dominant_lu) %>%
    mutate(max_hits = max(hits)) %>% 
    ungroup() %>%
    filter(hits == max_hits)
  
  length(unique(lu_worst$endPoint[lu_worst$dominant_lu %in% 'AgMix']))
  length(unique(lu_worst$endPoint[lu_worst$dominant_lu %in% 'Crops']))
  length(unique(lu_worst$endPoint[lu_worst$dominant_lu %in% 'Urban']))
  
  
    
  return(filt_top_ends)

}

