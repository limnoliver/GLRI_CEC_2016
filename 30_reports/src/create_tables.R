library(openxlsx)
library(toxEval)
library(dplyr)
library(tidyr)

get_siteTables <- function(chemicalSummary, chem_info, site_info, AOP, file_out){
  browser()
  chemicalSummary <- chemicalSummary %>%
    left_join(select(end_point_info, 
                     endPoint=assay_component_endpoint_name,
                     subFamily=intended_target_family_sub,
                     gene_symbol=intended_target_gene_symbol), by="endPoint") %>%
    left_join(select(chem_info, CAS, `Chemical Name`), by="CAS")
  
  chem_info <- chem_info %>%
    left_join(distinct(select(chemicalSummary, CAS, `Chemical Name`)), by=c("CAS", "Chemical Name")) %>%
    arrange(Class, `Chemical Name`)
  
  tableData <- chemicalSummary %>%
    rename(Chemical=`Chemical Name`,
           Family=Bio_category) %>%
    group_by(site, endPoint, Family, subFamily, gene_symbol, Chemical) 
  
  max_Samples <- tableData %>%
    summarize(sumEAR = sum(EAR)) %>% #Sum per date
    slice(which.max(sumEAR)) %>% # Gets max per date
    filter(sumEAR > 0) %>%
    data.frame() %>%
    spread(Chemical, sumEAR) %>%
    arrange(site, Family, subFamily, gene_symbol) %>%
    select(site, Family, subFamily, gene_symbol,endPoint, everything()) %>%
    mutate(maxSample = rowSums(.[-1:-5], na.rm = TRUE)) %>%
    select(site, Family, subFamily, gene_symbol, endPoint, maxSample)
  
  tableData <- tableData %>%
    summarize(maxEAR = max(EAR)) %>%
    filter(maxEAR > 0) %>%
    data.frame() %>%
    spread(Chemical, maxEAR) %>%
    arrange(site, Family, subFamily, gene_symbol) %>%
    select(site, Family, subFamily, gene_symbol,endPoint, everything()) %>%
    left_join(AOP,by="gene_symbol") %>%
    left_join(max_Samples, by=c("site", "Family", "subFamily", "gene_symbol","endPoint")) %>%
    select(site, Family, subFamily, gene_symbol, endPoint, AOP, maxSample, everything()) 
  
  list_tables <- list()
  
  for(i in 1:nrow(site_info)){
    
    site <- site_info$SiteID[i]
    site_name <- site_info$`Short Name`[i]
    tableData_site <- tableData[tableData$site == site,]
    tableData_site <- Filter(function(x)!all(is.na(x)), tableData_site)
    
    if(nrow(tableData_site) > 0){
      tableData2 <- select(tableData_site, -endPoint, -Family, -subFamily, -gene_symbol, -AOP, -site, -maxSample)
      tableData_site$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))
      orderedCols <- chem_info$`Chemical Name`[chem_info$`Chemical Name` %in% names(tableData_site)]
      
      tableData_site <- tableData_site[,c("Family", "subFamily", "gene_symbol", "endPoint","AOP", "maxSample", "nChems", orderedCols)] 
      
      list_tables[[site_name]] <- tableData_site
    }
    
  }
  
  write.xlsx(list_tables, file = file_out, append=TRUE)
  
}