library(readxl)
library(dplyr)

get_schedules <- function(schedule.numbers){
  
  info <- data.frame()
  
  for(schedule in schedule.numbers){

    url.info <- paste0("http://nwqlqc/servlets_u/AnalyticalServicesGuide?srchStr=",schedule,"&srchType=sched&mCrit=exact&oFmt=xl")
    temp.path <- tempdir()
    temp.file <- file.path(temp.path,"AnalyticalServicesGuide.xls")
    download.file(url.info, destfile = temp.file, mode="wb")
    schedule.data <- read_excel(temp.file)
    info <- bind_rows(info, schedule.data)
  }

  return(info)

}