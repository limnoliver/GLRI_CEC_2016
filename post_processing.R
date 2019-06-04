library(remake)
library(toxEval)

install_missing_packages()

make()
diagram()

tracking <- make("tracking")

chem_data <- make("chem_data")
chem_info <- make("chem_info")
chem_site <- make("site_info")

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = "intended_target_family_sub")

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info)
chemPlot <- plot_tox_boxplots(chemicalSummary, 
                              category = "Chemical",
                              plot_ND = FALSE)
chemPlot

explore_endpoints()
