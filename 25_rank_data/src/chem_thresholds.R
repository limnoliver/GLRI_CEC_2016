wq <- make('graph_data_wq')
tox <- make('graph_data_tox')

names(wq)[4] <- "bench"
names(tox)[4] <- "maxEAR"

compare <- full_join(wq, tox) %>%
  group_by(chnm, Class) %>%
  summarize(median_bench = median(bench, na.rm = T),
            median_maxEAR = median(maxEAR, na.rm = T))

cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
              "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e", "#800000", "#808000")

# calculate slopes fo lines where we have enough info to do so -
# that is, for herbicides and insecticides

herbicides <- filter(compare, Class %in% "Insecticide")
herb_lm <- lm(log10(herbicides$median_bench) ~ log10(herbicides$median_maxEAR))
herb_slope <- herb_lm$coefficients[2]

overall_lm <- lm(log10(compare$median_bench) ~ log10(compare$median_maxEAR))
ggplot(compare, aes(x = median_maxEAR, y = median_bench)) +
  geom_point(aes(color = Class)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(data = subset(compare, Class %in% c("Herbicide", "Insecticide")), 
              aes(group = Class, color = Class), method = 'lm', se = FALSE, alpha = 0.5) +
  scale_color_manual(values = cbValues[c(5,6,7,3,1,2,4)]) +
  theme_bw() +
  labs(x = "Median of max EAR per site", y = "Median of max TQ per site") +
  geom_abline(slope = overall_lm$coefficients[2], overall_lm$coefficients[1], alpha = 0.5) +
  geom_hline(yintercept = 0.1, color = 'black', linetype = 2, size = 1.5, alpha = 0.5)

ggsave()
c('orange', 'yellow', 'brown', 
  'green', 'red', 'blue', 'purple')
plot(compare$median_bench ~ compare$median_maxEAR,
     )
