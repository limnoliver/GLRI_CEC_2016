wq <- make('graph_data_wq') %>%
  filter(meanEAR > 0)
tox <- make('graph_data_tox') %>%
  filter(meanEAR > 0)

names(wq)[4] <- "bench"
names(tox)[3] <- "maxEAR"

compare <- full_join(wq, tox) %>%
  group_by(chnm, Class) %>%
  summarize(median_bench = median(bench, na.rm = T),
            median_maxEAR = median(maxEAR, na.rm = T))

ggplot(compare, aes(y = maxEAR, x = compare$chnm))+
  geom_boxplot() +
  scale_y_log10()

plot(compare$maxEAR ~ compare$bench, xlim = c(0,1), ylim = c(0,1))
cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
              "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e", "#800000", "#808000")

# calculate slopes fo lines where we have enough info to do so -
# that is, for herbicides and insecticides

herbicides <- filter(compare, Class %in% "Herbicide")
herb_lm <- lm(log10(herbicides$median_bench) ~ log10(herbicides$median_maxEAR))
herb_slope <- herb_lm$coefficients[2]

plot(log10(compare$median_maxEAR) ~ log10(compare$median_bench))
abline(lm(log10(compare$median_maxEAR) ~ log10(compare$median_bench)))
overall_lm <- lm(log10(compare$median_maxEAR) ~ log10(compare$median_bench))

p <- ggplot(compare, aes(x = median_bench, y = median_maxEAR)) +
  geom_point(aes(color = Class)) +
  scale_x_log10() +
  scale_y_log10() +
  #annotation_logticks() +
  #geom_smooth(data = subset(compare, Class %in% c("Herbicide", "Insecticide")), 
  #            aes(group = Class, color = Class), method = 'lm', se = FALSE, alpha = 0.5) +
  geom_smooth(method = 'lm', color = 'darkgray')+
  #geom_abline(slope = 1, intercept = 0, alpha = 0.8, color = 'darkgray') +
  scale_color_manual(values = cbValues[c(5,6,7,3,1,2,4)]) +
  theme_bw() +
  labs(x = "Median of max TQ per site", y = "Median of max EAR per site") +
  #geom_abline(slope = overall_lm$coefficients[2], overall_lm$coefficients[1], alpha = 0.5) +
  geom_vline(xintercept = 0.1, color = 'black', linetype = 2, size = 1, alpha = 0.5)

ggsave('figure/bench_vs_EAR.png', p, height = 3.5, width = 6)

