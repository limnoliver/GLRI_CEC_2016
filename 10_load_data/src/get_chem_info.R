get_chem_info <- function(infile) {
  info <- read.csv(infile, nrows = 240, stringsAsFactors = FALSE)
  
  info_c <- info %>%
    select(compound = Pesticide.compound, 
           pCode = USGS.Parameter.code, 
           CAS = CASRN, 
           parent_pesticide = Pesticide.group,
           MlWt = Molecular.weight) %>%
    mutate(pCode = as.character(pCode))
  
  usgs_pcodes_cas <- dataRetrieval::parameterCdFile %>%
    select(pCode = parameter_cd, CAS_2 = casrn) %>%
    filter(pCode %in% info_c$pCode)
  
  info_c <- left_join(info_c, usgs_pcodes_cas) %>%
    select(-CAS) %>%
    rename(CAS = CAS_2)
  
  # change incorrect CAS number
  # Hexazinone TP F	68574 should have cas == 56611-55-3
  
  info_c$CAS[info_c$pCode %in% '68574'] <- '56611-55-3'
  info_c$parent_pesticide[info_c$parent_pesticide %in% 'Bifenthrin; lambda-cyhalothrin; Tefluthrin'] <- 'Bifenthrin'

  # fix space after Azinphos-methyl
  info_c$parent_pesticide[grepl('Azinphos-methyl', info_c$parent_pesticide)] <- 'Azinphos-methyl'
  info_c$parent_pesticide[grepl('thiobencarb', info_c$parent_pesticide)] <- 'Thiobencarb'
  
  
  return(info_c)
  
  }