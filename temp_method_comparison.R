# compare data
# compare glyphosate and imidacloprid methods
library(dplyr)
library(tidyr)
library(ggplot2)
# glyphosate
head(glyphosate_clean)
gly_pcodes <- unique(glyphosate_clean$pCode)
parameterCdFile[parameterCdFile$parameter_cd %in% gly_pcodes, ]

# 62722 is other method for glyphosate
# 99960 is immunoassay for glyphosate
# 62649 is aminomethylphosphonic acid, a degradate of glyphosate

test <- filter(NWIS, pCode %in% '99960')

glyph_summary <- group_by(glyphosate_clean, SiteID, sample_dt) %>%
  summarize(n())

glyph_compar <- select(glyphosate_clean, SiteID, sample_dt, pCode, value, remark_cd) %>%
  #mutate(value = ifelse(remark_cd %in% "<", 0.5*value, value)) %>%
  group_by(SiteID, sample_dt, pCode) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  spread(key = pCode, value = value)

glyph_comments <- select(glyphosate_clean, SiteID, sample_dt, pCode, value, remark_cd) %>%
  #mutate(value = ifelse(remark_cd %in% "<", 0.5*value, value)) %>%
  group_by(SiteID, sample_dt, pCode) %>%
  summarize(remark_cd = ifelse(grep))

immuno_mdl <- filter(glyphosate_clean, pCode %in% '99960' & remark_cd %in% '<') %>%
  select(value) %>%
  distinct()

full_mdl <- filter(glyphosate_clean, pCode %in% '62722' & remark_cd %in% '<') %>%
  select(value) %>%
  distinct()
  

ggplot(glyph_compar, aes(x = `99960`, y = `62722`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(trans = 'log10', limits = c(0.01,1.5)) +
  scale_y_continuous(trans = 'log10', limits = c(0.01, 1.5)) +
  theme_bw()

ggplot(glyph_compar, aes(x = `62722`, y = `62649`)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_abline(slope = 2, intercept = 0, col = 'red') + 
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10')

# compare imidacloprid methods
imidacloprid_p <- filter(pesticides_clean, pCode == '68426')
imidacloprid_n <- filter(neonics_clean, pCode == '68426')  

imidacloprid <- bind_rows(imidacloprid_p, imidacloprid_n) %>%
  select(SiteID, sample_dt, value, source) %>%
  spread(key = source, value = value)

imidacloprid_c <- bind_rows(imidacloprid_p, imidacloprid_n) %>%
  select(SiteID, sample_dt, remark_cd, source) %>%
  spread(key = source, value = remark_cd) %>%
  rename(neonic_rmk = neonic, pesticides_s2437_rmk = pesticides_s2437)

imidacloprid <- left_join(imidacloprid, imidacloprid_c) %>%
  mutate(neonic = ifelse(neonic_rmk %in% '<', 0.5*neonic, neonic),
         pesticides_s2437 = ifelse(pesticides_s2437_rmk %in% '<', 0.5*pesticides_s2437, pesticides_s2437))

imid_rpd <- filter(imidacloprid, !is.na(neonic) & !is.na(pesticides_s2437)) %>%
  filter(!(neonic_rmk %in% "<") & !(pesticides_s2437_rmk %in% "<")) %>%
  mutate(rpd = 100*((neonic-pesticides_s2437)/((neonic + pesticides_s2437)/2))) %>%
  mutate(rpd_direction = ifelse(rpd <= 0, 'Pesticide > Neonic', 'Neonic > Pesticide'))

ggplot(imidacloprid, aes(x = pesticides_s2437, y = neonic)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'lm') + 
  theme_bw()

ggplot(imid_rpd, aes(x = pesticides_s2437, y = neonic)) +
  geom_point(alpha = 0.5, aes(size = abs(rpd), color = rpd_direction)) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'lm')

mdl_n <- filter(imidacloprid, neonic_rmk %in% '<') %>%
  select(neonic) %>%
  distinct()

mdl_p <- filter(imidacloprid, pesticides_s2437_rmk %in% '<') %>%
  select(pesticides_s2437) %>%
  distinct()


