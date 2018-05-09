library(dplyr)

get_AOPs <- function(file_name){
  AOP_crosswalk <- read.csv(file_name, stringsAsFactors = FALSE)
  AOP <- AOP_crosswalk %>%
    select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
    distinct()
  
  return(AOP)
}