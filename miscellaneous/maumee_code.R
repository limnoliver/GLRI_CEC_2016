# maumee data

chem_ear <- make('chemicalSummary_deg_meto') %>%
  filter(shortName %in% 'Maumee') %>%
  mutate(class2 = gsub('Deg - ', '', Class)) %>%
  group_by(parent_pesticide, date, class2) %>%
  summarize(parent_sum = sum(EAR)) %>%
  group_by(parent_pesticide, class2) %>%
  summarize(hits = length(which(parent_sum > 0.001)))

chem_bench <- make('chemicalSummary_bench_deg_meto') %>%
  filter(shortName %in% 'Maumee') %>%
  mutate(class2 = gsub('Deg - ', '', Class)) %>%
  group_by(CAS, parent_pesticide, date, class2) %>%
  summarize(max_tq = max(EAR)) %>%
  group_by(parent_pesticide, date, class2) %>%
  summarize(parent_sum = sum(max_tq)) %>%
  group_by(parent_pesticide, class2) %>%
  summarize(hits = length(which(parent_sum > 0.1))) %>%
  arrange(class2, hits)

# exclusions luke already went through:
done_exc <- readxl::read_xlsx('data/Exclusions_toxEval32.xlsx', sheet = 'Notes')

# detected chems with EAR values
detected_chems <- make('chemicalSummary') %>%
  select(CAS, chnm) %>%
  distinct()

chems_bench <- make('chemicalSummary_bench')%>%
  select(CAS, chnm) %>%
  distinct()

chem_endpoints <- make('chemicalSummary') %>%
  select(CAS, chnm, endPoint) %>%
  distinct()

acc <- toxEval::ToxCast_ACC
needs_exc <- filter(detected_chems, !CAS %in% done_exc$CAS)

met <- filter(chem_endpoints, CAS %in% '51218-45-2')
met_acc <- filter(acc, CAS %in% '51218-45-2')

ory <- filter(chem_endpoints, CAS %in% '19044-88-3')
ory_acc <- filter(acc, CAS %in% '19044-88-3')

pros <- filter(make('chemicalSummary'), CAS %in% '94125-34-5')
pros_acc <- filter(acc, CAS %in% '94125-34-5')

pip <-  filter(chem_endpoints, CAS %in% '51-03-6')
pip_acc <- filter(acc, CAS %in% '51-03-6')

View(filter(met_acc, endPoint %in% met$endPoint))
