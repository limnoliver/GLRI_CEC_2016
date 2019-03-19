get_chem_info <- function(infile) {
  info <- read.csv(infile, nrows = 240, stringsAsFactors = FALSE)
  
  info_c <- info %>%
    select(compound = Pesticide.compound, 
           pCode = USGS.Parameter.code, 
           CAS = CASRN, 
           parent_pesticide = Pesticide.group,
           MlWt = Molecular.weight)
  
  return(info_c)
  
  }