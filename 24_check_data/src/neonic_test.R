# check neonic tox

summarize_neonic_tox <- function(tox_dat) {
 
  pcode_info <- make('pCodeInfo', remake_file= '10_load_data.yml')
   neonics <- make('neonic', remake_file = '10_load_data.yml') %>%
    select(pCode) %>%
    distinct() %>%
    left_join(pcode_info, by = c('pCode'= 'parameter_cd'))
  
  EAR_sums <- chemicalSummary %>% 
    group_by(site, date, chnm, CAS, Class) %>% 
    summarize(sumEAR = sum(EAR)) %>%
    filter(CAS %in% neonics$casrn) %>%
    filter(sumEAR > 0.001)
  
  conc_dat <- make('chemicalSummary_conc', remake_file = '20_merge_data.yml') %>%
    filter(CAS %in% neonics$casrn)
  
  al_bench <- make('graph_data_wq', remake_file = '40_pesticide_figs.yml') %>%
    filter(chnm %in% as.character(unique(conc_dat$chnm))) %>%
    filter(meanEAR > 0.1)
  
  hbs_bench <- make('graph_data_hbs', remake_file = '40_pesticide_figs.yml') %>%
    filter(chnm %in% as.character(unique(conc_dat$chnm))) %>%
    filter(meanEAR > 0.1)

}