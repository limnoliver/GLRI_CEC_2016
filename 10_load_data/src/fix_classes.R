fix_classes <- function(raw_classes) {
  
  #fipronil amide currently classified as 'other'
  #when it should be insecticide degradate
  raw_classes$Class[raw_classes$CAS == '205650-69-7'] <- 'Deg - Insecticide'
  raw_classes$Class[raw_classes$CAS == '2327-02-8'] <- 'Deg - Herbicide'
  raw_classes$Class[raw_classes$CAS == '35045-02-4'] <- 'Deg - Herbicide'
  
  #Desulfinylfipronil amide
  raw_classes$Class[raw_classes$CAS == '1115248-09-3'] <- 'Deg - Insecticide'
  
  #EPTC and degradates as "herbicides"
  raw_classes$Class[raw_classes$CAS == '759-94-4'] <- 'Herbicide'
  raw_classes$Class[raw_classes$CAS == '65109-69-5'] <- 'Deg - Herbicide'
  
  # Tebuthiuron degradate should be herbicide deg
  raw_classes$Class[raw_classes$CAS == '139888-73-6'] <- 'Deg - Herbicide'
  
  # 1H-1,2,4-Triazole as degradate of Propiconizole
  raw_classes$Class[raw_classes$CAS == '288-88-0'] <- 'Deg - Fungicide'
  
  # Dichlorovos can be a product and degradate. Calling it a degradate. 
  raw_classes$Class[raw_classes$CAS == '62-73-7'] <- 'Deg - Insecticide'
  
  # triallate metabolite 2,3,3-trichloro-prop-2-en-sulfonic acid
  raw_classes$Class[raw_classes$CAS == '65600-61-5'] <- 'Deg - Herbicide'
  
  raw_classes <- add_row(raw_classes, CAS = '56611-55-3', Class = 'Deg - Herbicide')
  
  return(raw_classes)
}