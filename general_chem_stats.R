# compare site rankings based on benchmarks vs EAR

site_rankings <- make('site_rankings')
library(dplyr)
library(tidyr)
site_rank_wide <- select(site_rankings, site, date, median_rank, variable) %>%
  spread(key = variable, value = median_rank)

plot(site_rank_wide$max_bench ~ site_rank_wide$max_sumEAR)
plot(site_rank_wide$mean_conc ~ site_rank_wide$n_detected)

abline(0,1)

plot(log(site_rank_wide_vals$mean_bench) ~ log(site_rank_wide_vals$mean_sumEAR))


head(site_rank_wide)

# which chems in toxCast have benchmarks?
conc <- select(chemicalSummary_conc, CAS, Class) %>% distinct()
bench <- select(chemicalSummary_bench, CAS, Class) %>% distinct()
toxcast <- select(chemicalSummary, CAS, Class) %>% distinct()

length(which(toxcast %in% bench))
length(which(bench %in% toxcast))

head(chemicalSummary)
