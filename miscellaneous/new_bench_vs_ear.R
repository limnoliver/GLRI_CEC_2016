ears <- make('chemicalSummary') %>%
  filter(EAR > 0) %>%
  group_by(CAS, chnm, site, date, Class) %>%
  summarize(sumEAR = sum(EAR))

chems_observed <- make('chemicalSummary_conc')

acc <- toxEval::ToxCast_ACC %>%
  filter(CAS %in% chems_observed$CAS) %>% 
  group_by()

conc/acc = ear



bench <- make('chemicalSummary_bench') %>%
  filter(EAR > 0) %>%
  group_by(CAS, chnm, site, date, Class) %>%
  summarize(maxBench = max(EAR))

compare <- full_join(ears, bench) %>%
  mutate(ratio = sumEAR/maxBench) %>%
  group_by(CAS, chnm, Class) %>%
  summarize(median_maxBench = median(maxBench),
            median_sumEAR = median(sumEAR))


plot(compare$median_sumEAR ~ compare$median_maxBench, xlim = c(0,1), ylim = c(0,1))
cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
              "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e", "#800000", "#808000")

# calculate slopes fo lines where we have enough info to do so -
# that is, for herbicides and insecticides

herbicides <- filter(compare, Class %in% "Herbicide") %>%
  filter(!is.na(median_maxBench)) %>% filter(!is.na(median_sumEAR))
herb_lm <- lm(log10(herbicides$median_maxBench) ~ log10(herbicides$median_sumEAR))
herb_slope <- herb_lm$coefficients[2]

plot(log10(compare$median_maxEAR) ~ log10(compare$median_bench))
abline(lm(log10(compare$median_maxEAR) ~ log10(compare$median_bench)))
overall_lm <- lm(log10(compare$median_maxEAR) ~ log10(compare$median_bench))

compare <- compare %>%
  mutate(type = ifelse(grepl('deg', Class, ignore.case = TRUE), 'degradate', 'parent')) %>%
  mutate(Class2 = case_when(
    grepl('Fungicide', Class) ~ 'Fungicide',
    grepl('Herbicide', Class) ~ 'Herbicide',
    grepl('Insecticide', Class) ~ 'Insecticide',
    TRUE ~ 'Other'
  )) %>%
  filter(!is.na(median_maxBench)) %>%
  filter(!is.na(median_sumEAR))

p <- ggplot(compare, aes(x = median_maxBench, y = median_sumEAR)) +
  geom_point(aes(color = Class2, shape = type)) +
  geom_hline(yintercept = 0.001, color = 'red') +
  
  scale_x_log10() +
  scale_y_log10() +
  #annotation_logticks() +
  #geom_smooth(data = subset(compare, Class %in% c("Herbicide", "Insecticide")), 
  #            aes(group = Class, color = Class), method = 'lm', se = FALSE, alpha = 0.5) +
  geom_smooth(method = 'lm', color = 'darkgray')+
  #geom_abline(slope = 1, intercept = 0, alpha = 0.8, color = 'darkgray') +
  #scale_color_manual(values = cbValues[c(1,2,3,1,2,3,4)]) +
  scale_shape_manual(values = c(1, 16))+
  theme_bw() +
  labs(x = "Median of max TQ per site", y = "Median of sum EAR per site", color = 'Class', shape = 'Type') +
  #geom_abline(slope = overall_lm$coefficients[2], overall_lm$coefficients[1], alpha = 0.5) +
  geom_vline(xintercept = 0.1, color = 'black', linetype = 2, size = 1, alpha = 0.5)

ggsave('figures/ms_supplement_figures/bench_vs_EAR.png', p, height = 4.7, width = 6)
