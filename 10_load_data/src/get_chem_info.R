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
  
  # filter a few compounds out that have repeated CAS numbers
  #pCode 68566 (Demethyl hexazinone B) and Hexazinone TP F	68574 have same CAS
  
  info_c <- filter(info_c, !pCode == '68574')
  
  # fix missing parent/compounds
  
  return(info_c)
  
  }