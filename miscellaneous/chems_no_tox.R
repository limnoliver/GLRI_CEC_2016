# which compounds aren't in toxcast or benchmarks

conc <- make('chemicalSummary_conc')
bench <- make('chemicalSummary_bench')
bench_deg <- make('chemicalSummary_bench_deg_meto')
ear <- make('chemicalSummary')
info <- make('chem_master')

# test metolachlor OA 152019-73-3

filter(bench, CAS %in% '152019-73-3')
filter(bench_deg, CAS %in% '152019-73-3')
filter(conc, CAS %in% '152019-73-3')

filter(ear, CAS %in% '152019-73-3')

length(unique(conc$CAS))
length(unique(bench$CAS))
length(unique(ear$CAS))

all_cas <- unique(conc$CAS)[!unique(conc$CAS) %in% unique(bench$CAS, ear$CAS)]

out <- filter(info, CAS %in% all_cas) %>%
  mutate(CAS = ifelse(grepl('fake', CAS), '', CAS))

write.csv(out, 'chems_no_bench_toxcast.csv', row.names = FALSE)
