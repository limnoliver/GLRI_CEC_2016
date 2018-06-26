clean_pesticides <- function(dat) {

  dat$source <- 'pesticides_s2437'
  return(dat)
}

clean_neonics <- function(dat) {
  
  dat$source <- 'neonic'
  return(dat)
}

clean_glyphosate <- function(dat) {
  dat$source <- 'glyphosate'
  return(dat)
}