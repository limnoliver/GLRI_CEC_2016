library(dplyr)
library(USGSHydroTools)

merged_neonic_flow <- function(neonic, flow, sites){

  flow$Q <- ifelse(!is.na(flow$X_00060_00000),flow$X_00060_00000,flow$X_00060_00003)
  
  listQ <- list()
  
  for(site in sites$USGS.station.number){
    
    if(site == "04157005") site <- "04157000"
    
    sub_data <- filter(neonic, USGS.station.number == site)
    dfQ2 <- filter(flow, site_no == site)
    
    if(nrow(dfQ2) > 0){
      subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = sub_data, starttime = "pdate",
                             times=6,units="hrs",stats.return = c("mean","max","min","difference"),
                             subdatesvar = "",out.varname = "Q")
      subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                             times=0,units="hrs",stats.return = c("nearprev"),
                             subdatesvar = "",out.varname = "Q")
      subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                             times=1,units="hrs",stats.return = c("mean"),
                             subdatesvar = "",out.varname = "Q")
      listQ[[site]] <- subdfNeonic
    }
    
  }

  dfNeonic2 <- listQ[[1]]

  for(i in 2:length(listQ)){
    dfNeonic2 <- bind_rows(dfNeonic2,listQ[[i]])
  }
  
  return(dfNeonic2)
}

merged_NWIS <- function(tracking, NWIS, neonic_flow){
  
  browser()
  
}
