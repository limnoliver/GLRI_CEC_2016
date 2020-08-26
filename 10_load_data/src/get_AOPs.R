library(dplyr)

get_AOPs <- function(file_name){
  AOP_crosswalk <- read.csv(file_name, stringsAsFactors = FALSE)
  AOP <- AOP_crosswalk %>%
    select(endPointID = `Assay.Endpoint.ID`,
           endPoint=Component.Endpoint.Name, 
           AOPID=AOP..,
           AOPTitle = `AOP.Title`,
           KeyID = `KE.`,
           KeyName = `Key.Event.Name`,
           KeyType = `KeyEvent.Type`) %>%
    distinct()
  
  return(AOP)
}