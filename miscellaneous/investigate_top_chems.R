# investigate diazinon and alachlor
ear <- make('chemicalSummary_deg_meto')
bench <- make('chemicalSummary_bench_deg_meto')

# diazinon
diaz_ear <- filter(ear, parent_pesticide %in% 'Diazinon') %>% 
  group_by(chnm, site, date) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(chnm) %>%
  summarize(detects = n(),
            hits = length(which(sumEAR > 0.001)))

diaz_bench <- filter(bench, parent_pesticide %in% 'Diazinon') %>% 
  group_by(chnm, site, date) %>%
  summarize(maxBench = max(EAR)) %>%
  group_by(chnm) %>%
  summarize(hits = length(which(maxBench > 0.1)))

