library(dplyr)
library(dataRetrieval)

get_flow <- function(tracking){
  
  sitesIDS <- zeroPad(unique(tracking$SiteID),padTo = 8)
  begin <- "2015-10-01"
  end <- "2016-10-01"
  timeSteps <- as.numeric((difftime(as.POSIXct(end),as.POSIXct(begin),units="days")+2)*96)
  datesVector <- as.POSIXct(begin,tz="UTC") + 5*60*60 + seq(from = 1, to = timeSteps)*15*60
  dfDates <- data.frame(dateTime=datesVector)
  
  dfQ <- readNWISuv(sitesIDS, "00060", begin, end)
  dfDMQ <- readNWISdv(sitesIDS, "00060", begin, end)
  
  dfQ1 <- right_join(dfQ,dfDates,by = c("dateTime"))
  dfQ1$Date <- as.Date(dfQ1$dateTime)
  dfQ2 <- left_join(dfQ1,dfDMQ,by=c("Date","agency_cd","site_no"))
  dfQ2$Q <- ifelse(!is.na(dfQ2$X_00060_00000),dfQ2$X_00060_00000,dfQ2$X_00060_00003)
  
  return(dfQ2)
  

}