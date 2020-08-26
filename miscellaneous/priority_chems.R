# top chems
chem_master <- make('chem_master')
# all measurements
chem_count <- make('chem_data_complete') %>%
  group_by(pCode) %>%
  summarize(n_measured = n()) %>%
  left_join(select(chem_master, pCode, CAS, Class, parent_pesticide))

chem_detect <- make('chemicalSummary_conc') %>%
  group_by(chnm, CAS) %>%
  summarize(n_detect = n(), n_sites = length(unique(site))) %>%
  left_join(select(chem_count, CAS, n_measured, Class, parent_pesticide)) %>%
  mutate(detect_rate = n_detect/n_measured)

chem_detect_top <- chem_detect %>%
  filter(detect_rate >= 0.25) %>%
  ungroup() %>% pull(CAS)


chem_ear_hit <- make('chemicalSummary') %>%
  group_by(site, date, chnm, CAS) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(chnm, CAS) %>%
  summarize(n_ear_hit = length(which(sumEAR > 0.001))) %>%
  left_join(chem_count) %>%
  mutate(ear_hit_rate = n_ear_hit/n_measured)

chem_ear_top <- chem_ear_hit %>%
  filter(ear_hit_rate > 0)%>%
  ungroup() %>% pull(CAS)

chem_bench_hit <- make('chemicalSummary_bench') %>%
  group_by(site, date, chnm, CAS) %>%
  summarize(sumEAR = max(EAR)) %>%
  group_by(chnm, CAS) %>%
  summarize(n_bench_hit = length(which(sumEAR > 0.1))) %>%
  left_join(chem_count) %>%
  mutate(bench_hit_rate = n_bench_hit/n_measured)

chem_bench_top <- chem_bench_hit %>%
  filter(bench_hit_rate > 0) %>%
  ungroup() %>% pull(CAS)

# now calculate stats for priority chems
top <-  unique(c(chem_detect_top, chem_ear_top, chem_bench_top))

sample_chem_ear <- make('chemicalSummary') %>%
  filter(CAS %in% top) %>%
  group_by(site, date, chnm, CAS) %>%
  summarize(sumEAR = sum(EAR)) %>% ungroup()

sample_chem_bench <- make('chemicalSummary_bench') %>%
  filter(CAS %in% top) %>%
  group_by(site, date, chnm, CAS) %>%
  summarize(maxBench = max(EAR)) %>% ungroup()

sample_detect <- make('chemicalSummary_conc') %>%
  filter(CAS %in% top) %>%
  group_by(site, date, chnm, CAS) %>%
  summarize(detected = TRUE) %>% ungroup()

chem_count <- make('chem_data_complete') %>%
  left_join(select(chem_master, pCode, CAS)) %>%
  filter(CAS %in% top) %>%
  group_by(SiteID, `Sample Date`, CAS) %>%
  summarize(measured = TRUE) %>%
  select(site=SiteID, date=`Sample Date`, CAS, measured) %>% ungroup() %>%
  mutate(date = as.POSIXct(date))

all_samples <- make('chem_data_complete') %>%
  select(site=SiteID, date=`Sample Date`) %>%
  mutate(date = as.POSIXct(date)) %>%
  distinct(site, date) %>% 
  mutate(sample_id = 1:nrow(.))

all_possible <- expand.grid(CAS = top, sample_id = all_samples$sample_id) %>%
  left_join(all_samples)

  
samples <- all_possible %>%
  left_join(chem_count) %>%
  left_join(sample_detect) %>%
  left_join(sample_chem_ear) %>%
  left_join(sample_chem_bench) %>%
  mutate(measured = ifelse(is.na(measured), FALSE, TRUE),
         detected = ifelse(is.na(detected), FALSE, TRUE))

sample_names <- filter(chem_master, CAS %in% top) %>%
  select(CAS, compound, `Chemical Name`, Class) %>%
  mutate(compound = ifelse(is.na(compound), `Chemical Name`, compound)) %>%
  mutate(compound = ifelse(CAS == '19988-24-0', 'Deethylhydroxyatrazine', compound)) %>%
  mutate(compound = ifelse(CAS == '2163-68-0', 'Hydroxyatrazine', compound)) %>%
  mutate(compound = ifelse(CAS == '3567-62-2', 'Monomethyldiuron', compound)) %>%
  mutate(compound = ifelse(grepl('Deg', Class), paste0('*', compound), compound)) %>%
  mutate(Class = gsub('Deg - ', '', Class))
  
  
# change class to parent class, also put asterisk next to degs

chem_summary <- samples %>%
  mutate(category = case_when(
    sumEAR > 0.01 | maxBench > 1 ~ 'g_hit_plus',
    (sumEAR > 0.001 & sumEAR <= 0.01) | (maxBench > 0.1 & maxBench <= 1) ~ 'f_hit',
    (sumEAR > 0.0001 & sumEAR <= 0.001) | (maxBench > 0.01 & maxBench <= .1) ~ 'e_hit_minus',
    sumEAR <= 0.0001 | maxBench <= .01 ~ 'd_detect',
    is.na(sumEAR) & is.na(maxBench) & detected == TRUE ~ 'c_detect_no_info',
    measured & !detected  ~ 'b_not_detected',
    !measured ~ 'a_not_measured'
  )) %>%
  group_by(CAS, category) %>%
  summarize(n = n()) %>%
  left_join(select(sample_names, CAS, compound, Class)) 

chem_summary$Class <- factor(chem_summary$Class, levels = c('Herbicide', "Insecticide", 'Fungicide'))

# find order of chems within class
chem_order <- ungroup(chem_summary) %>%
  #filter(category %in% c('g_hit_plus', 'f_hit')) %>%
  group_by(compound, CAS) %>%
  summarize(value = sum(n[category %in% c('g_hit_plus', 'f_hit')]),
            value2 = sum(n[category %in% c('e_hit_minus', 'd_detect')]),
            value3 = sum(n[category %in% 'c_detect_no_info'])) %>%
  arrange(value, value2, value3)

chem_summary$compound <- factor(chem_summary$compound, levels = chem_order$compound)
  

#colors for stack bar
colors_EAR <- brewer.pal(n = 9, name = "YlOrRd")[c(2,4,7,9)]

#Color if using a grey to show unknowns
colors_EAR2 <- c('grey90', 'grey60', 'grey30', colors_EAR)

p <- ggplot(data=chem_summary, aes(x=compound, y=n, fill=category)) + 
  geom_bar(color = 'gray10', width=.8, size=.3, stat='identity') +  
  coord_flip() +
  facet_grid(rows = vars(Class), space="free", scales="free") +
  labs(x=NULL, y='Number of samples', fill = '') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'bottom', 
        strip.background = element_blank()) + 
  guides(fill = guide_legend(ncol = 2, byrow = FALSE, reverse = TRUE, label.hjust = 0)) +
  scale_fill_manual(values = colors_EAR2,
                    labels = c("Not measured",
                               "Not detected",
                               "Detected, no benchmark",
                               expression(paste("EAR < 10"^"-4", " and TQ < 10"^"-2")),
                               expression(paste("EAR > 10"^"-4", " or TQ > 10"^"-2")),
                               expression(paste("EAR > 10"^"-3", " or TQ > 10"^"-1")),
                               expression(paste("EAR > 10"^"-2", " or TQ > 1")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(expand=c(0,2))  
  #theme(legend.position = 'bottom', legend.text.align = 0) +
  #guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T, nrow=2)) +
  #theme(axis.text.x = element_text(vjust=0.5, hjust=1))
  

ggsave('figures/ms_figures/top_chems_stacked_barplot.png', p, height = 8, width = 4.5)
