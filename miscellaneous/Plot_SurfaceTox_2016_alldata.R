
#Plot toxeval results for Oliver et al. 2020
#Includes data from 16 rivers/streams in Great Lakes
#biweekly to monthly sampling for one year

#tox eval summaries, note these are missing chemicals not in toxcast
chemicalSummary_surface #Sam's data from surface
chemicalSummary_bench_surface # Sam's surface data using custom benchmarks

#This is the tox_list file that was used to generate the chemical summary
# It contains all chemicals
tox_list_surface$chem_data # Data prior to tox eval

#Check out the standard boxplots
plot_tox_boxplots(chemicalSummary_bench_surface, 
                  category = "Chemical", 
                  sum_logic = FALSE,
                  x_label = "Toxicity Quotient")

plot_tox_boxplots(chemicalSummary_surface, 
                  category = "Chemical", 
                  sum_logic = FALSE)

#colors for stack bar
colors_EAR <- brewer.pal(n = 9, name = "YlOrRd")[c(2,4,7,9)]

#Color if using a grey to show unknowns
colors_EAR2 <- c('grey90', colors_EAR)

#Figure out how many total chemicals were detected in each site
chnm_list <- tox_list_surface$chem_data %>%
  filter (Value > 0) %>%
  group_by (SiteID, CAS, pCode) %>%
  summarize(Value = max(Value, na.rm=T)) %>%
  group_by(SiteID) %>%
  summarize(n=n()) %>%
  left_join(tox_list_surface$chem_site) %>%
  rename(shortName = `Short Name`)

saveOutput = write.csv(chnm_list, file = file_out(file.path(path_to_data, "Data", "SurfaceWaterChemDetects_AllDates.csv")), row.names = F)


# #######################
# Summarize EAR by site
# #######################


chemicalSummary2_surface <- chemicalSummary_surface %>%
  # filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(EAR>0) %>%
  # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  # filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  dplyr::group_by(site, chnm, date) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>% #Take Max for each chnm, site, and date
  dplyr::group_by(site, chnm) %>%
  summarize(EAR=max(EAR, na.rm=T)) #Take max among dates


# Count number of chemicals per stream at each EAR level
site_freq_EAR0.01<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.01")

site_freq_EAR0.001<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.001 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.001")

site_freq_EAR0.0001<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR < 0.001, EAR >= 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveEAR0.0001")

site_freq_EARDetected<-chemicalSummary2_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR > 0 & EAR < 0.0001) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "EARDetected")

#Order of sites and names for plotting
site_table <- unique(chemicalSummary_surface[c('shortName', 'site')])
site_order_surface <- site_table$shortName[match(surface_ID_order, site_table$site)]

site_detection <- full_join(site_freq_EAR0.01, site_freq_EAR0.001) %>%
  full_join(site_freq_EAR0.0001) %>%
  full_join(site_freq_EARDetected) %>%
  left_join(unique(chemicalSummary_surface[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  left_join(chnm_list[c('shortName', 'n')]) %>%
  mutate(Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order_surface),
         UnknownTox = n - sum(AboveEAR0.01, AboveEAR0.001, AboveEAR0.0001, EARDetected, na.rm=T))

site_detection <- site_detection %>%
  arrange(Lake, shortName) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake, -n) %>%
  mutate(group = factor(group, levels=rev(c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "EARDetected", "UnknownTox"))))

site_detection$value[which(is.na(site_detection$value))] <- 0


#Plot stacked barplot by site 
sitebyEAR_surfaceall <- ggplot(data=site_detection, aes(x=shortName, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Lake, space="free", scales="free") +
  labs(x='Stream', y='Number of chemicals', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR2, labels = c("EAR Unknown", expression(paste("EAR < 10"^"-4")), expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")), expression(paste("EAR > 10"^"-2")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(expand=c(0,1)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
print(sitebyEAR_surfaceall)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_EARBySite_watersamples_allyear.png")), plot = sitebyEAR_surfaceall, height=4, width=6)



#Make similar table for TQ's

chemicalSummary2_bench_surface <- chemicalSummary_bench_surface %>%
  # filter(date>=date_filter[1], date<=date_filter[2]) %>%
  filter(EAR>0) %>%
  # filter(CAS %in% unique(tox_list_allpocis$chem_data$CAS)) %>%
  # filter(site %in% unique(tox_list_allpocis$chem_site$SiteID)) %>%
  dplyr::group_by(site, chnm, date) %>% 
  summarize(EAR = max(EAR, na.rm=T)) %>% #Take Max for each chnm, site, and date
  dplyr::group_by(site, chnm) %>%
  summarize(EAR=max(EAR, na.rm=T)) #Take max among dates


# Count number of chemicals per stream at each EAR level
site_freq_TQ1<-chemicalSummary2_bench_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ1")

site_freq_TQ0.1<-chemicalSummary2_bench_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.1 & EAR < 1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ0.1")

site_freq_TQ0.01<-chemicalSummary2_bench_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0.01 & EAR < 0.1) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "AboveTQ0.01")

site_freq_TQDetected<-chemicalSummary2_bench_surface %>%
  group_by(site, chnm) %>%
  dplyr::filter(EAR >= 0 & EAR < 0.01) %>%
  summarize(EAR_max = max(EAR, na.rm=T)) %>%
  tally(name = "TQDetected")


site_detection_TQ <- full_join(site_freq_TQ1, site_freq_TQ0.1) %>%
  full_join(site_freq_TQ0.01) %>%
  full_join(site_freq_TQDetected) %>%
  left_join(unique(chemicalSummary_surface[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  left_join(chnm_list[c('shortName', 'n')]) %>%
  mutate(Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order_surface),
         UnknownTox = n - sum(AboveTQ1, AboveTQ0.1, AboveTQ0.01, TQDetected, na.rm=T))



site_detection_TQ <- site_detection_TQ %>%
  arrange(Lake, shortName) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake, -n) %>%
  mutate(group = factor(group, levels=rev(c("AboveTQ1", "AboveTQ0.1", "AboveTQ0.01", "TQDetected", "UnknownTox"))))

site_detection_TQ$value[which(is.na(site_detection_TQ$value))] <- 0




#Plot barplot by site 
sitebyTQ_surfaceall <- ggplot(data=site_detection_TQ, aes(x=shortName, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Lake, space="free", scales="free") +
  labs(x='Stream', y='Number of chemicals', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR2, labels = c("Unknown TQ",expression(paste("TQ < 10"^"-2")), expression(paste("TQ > 10"^"-2")), expression(paste("TQ > 10"^"-1")), expression(paste("TQ > 1")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(expand=c(0,1)) + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
print(sitebyTQ_surfaceall)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_TQBySite_watersamples_allyear.png")), plot = sitebyTQ_surfaceall, height=4, width=6)



#combine EAR and TQ

EAR_thresholds <-data.frame(group=c("AboveEAR0.01", "AboveEAR0.001", "AboveEAR0.0001", "EARDetected"),
                            min=c(0.01, 0.001, 0.0001, 0), 
                            max=c(Inf, 0.01, 0.001, 0.0001))

TQ_thresholds <-data.frame(group=c("AboveTQ1", "AboveTQ0.1", "AboveTQ0.01", "TQDetected"),
                            min=c(1, 0.1, 0.01, 0), 
                            max=c(Inf, 1, 0.1, 0.01))

#Divide TQ by 100 and take max of TQ or EAR. 
chemicalSummary2_combine <- chemicalSummary2_bench_surface %>%
  rename(TQ = EAR) %>%
  full_join(chemicalSummary2_surface) %>%
  mutate(TQ_X2 = TQ/100,
    TQ_group = NA,
         EAR_group=NA) %>%
  rowwise() %>%
      mutate(Max = max(c(EAR, TQ_X2), na.rm=T))


# Count number of chemicals per stream at each EAR level
site_freq_Max0.01<-chemicalSummary2_combine %>%
  group_by(site, chnm) %>%
  dplyr::filter(Max >= 0.01) %>%
  summarize(EAR_max = max(Max, na.rm=T)) %>%
  tally(name = "AboveMax0.01")

site_freq_Max0.001<-chemicalSummary2_combine %>%
  group_by(site, chnm) %>%
  dplyr::filter(Max >= 0.001 & Max < 0.01) %>%
  summarize(EAR_max = max(Max, na.rm=T)) %>%
  tally(name = "AboveMax0.001")

site_freq_Max0.0001<-chemicalSummary2_combine %>%
  group_by(site, chnm) %>%
  dplyr::filter(Max >= 0.0001 & Max < 0.001) %>%
  summarize(EAR_max = max(Max, na.rm=T)) %>%
  tally(name = "AboveMax0.0001")

site_freq_MaxDetected<-chemicalSummary2_combine %>%
  group_by(site, chnm) %>%
  dplyr::filter(Max > 0 & Max < 0.0001) %>%
  summarize(EAR_max = max(Max, na.rm=T)) %>%
  tally(name = "MaxDetected")



site_detection_max <- full_join(site_freq_Max0.01, site_freq_Max0.001) %>%
  full_join(site_freq_Max0.0001) %>%
  full_join(site_freq_MaxDetected) %>%
  left_join(unique(chemicalSummary_surface[c("site", "shortName", "Lake")])) %>%
  rowwise() %>%
  left_join(chnm_list[c('shortName', 'n')]) %>%
  mutate(Lake = factor(Lake, c("Superior", "Michigan", "Huron", "Erie", "Ontario")),
         shortName = factor(shortName, site_order_surface),
         UnknownTox = n - sum(AboveMax0.01, AboveMax0.001, AboveMax0.0001, MaxDetected, na.rm=T))



site_detection_max <- site_detection_max %>%
  arrange(Lake, shortName) %>%
  tidyr::gather(key=group, value=value, -site, -shortName, -Lake, -n) %>%
  mutate(group = factor(group, levels=rev(c("AboveMax0.01", "AboveMax0.001", "AboveMax0.0001", "MaxDetected", "UnknownTox"))))

site_detection_max$value[which(is.na(site_detection_max$value))] <- 0





#Plot barplot by site 
sitebyMax_surfaceall <- ggplot(data=site_detection_max, aes(x=shortName, y=value, fill=group)) + 
  geom_bar(color = 'grey', width=.8, size=.1, stat='identity') +  
  # coord_flip() +
  facet_grid(.~Lake, space="free", scales="free") +
  labs(x='River', y='Number of chemicals', fill = 'EAR') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=.5),
        legend.title = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8)) + 
  scale_fill_manual(values = colors_EAR2,
                    labels = c("Unknown EAR and TQ",
                               expression(paste("EAR < 10"^"-4", " and TQ < 10"^"-2")),
                               expression(paste("EAR > 10"^"-4", " or TQ > 10"^"-2")),
                               expression(paste("EAR > 10"^"-3", " or TQ > 10"^"-1")),
                               expression(paste("EAR > 10"^"-2", " or TQ > 1")))) +
  # scale_fill_brewer(palette = "YlOrRd", labels = c("Detected", expression(paste("EAR > 10"^"-4")), expression(paste("EAR > 10"^"-3")))) +
  scale_y_continuous(expand=c(0,1)) + 
  theme(legend.position = 'bottom', legend.text.align = 0) +
  guides(fill = guide_legend(title.position='left', title.hjust=0.5, reverse=T, nrow=2)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
print(sitebyMax_surfaceall)

ggsave(file_out(file.path(path_to_data, "Figures/StackBox_EARorTQ_BySite_watersamples_allyear.png")), plot = sitebyMax_surfaceall, height=4, width=6)


