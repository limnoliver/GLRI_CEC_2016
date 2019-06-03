fix_classes <- function(raw_classes) {
  
  #fipronil amide currently classified as 'other'
  #when it should be insecticide degradate
  raw_classes$Class[raw_classes$CAS == '205650-69-7'] <- 'Deg - Insecticide'
  raw_classes$Class[raw_classes$CAS == '2327-02-8'] <- 'Deg - Herbicide'
  raw_classes$Class[raw_classes$CAS == '35045-02-4'] <- 'Deg - Herbicide'
  
  return(raw_classes)
}