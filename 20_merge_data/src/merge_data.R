library(dplyr)
library(USGSHydroTools)
library(tidyr)

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

merged_NWIS <- function(tracking, NWIS, neonic_flow, pCodeInfo){
  
  just_neonic_data <- neonic_flow %>%
    select(Sample, site = USGS.station.number, pdate, 
           Acetamiprid, 
           Clothianidin,
           Dinotefuran,
           Imidacloprid,
           Thiacloprid,
           Thiamethoxam) %>%
    gather(chemical, value, -Sample, -site, -pdate)
  
  just_neonic_remarks <- neonic_flow %>%
    select(Sample, site = USGS.station.number, pdate, 
           R_Acetamiprid, 
           R_Clothianidin,
           R_Dinotefuran,
           R_Imidacloprid,
           R_Thiacloprid,
           R_Thiamethoxam
    ) %>%
    gather(chemical_rk, remark_cd, -Sample, -site, -pdate) %>%
    mutate(chemical = gsub("R_","",chemical_rk)) %>%
    select(-chemical_rk)

  just_neonic <- left_join(just_neonic_data, 
                           just_neonic_remarks, by=c("Sample","site","pdate","chemical"))
  
  just_NWIS <- select(NWIS, site=SiteID, Sample = NWISRecordNumber, pdate, pCode, value) %>%
    left_join(select(pCodeInfo, pCode=parameter_cd, chemical=casrn), by="pCode") 
  
  nwis_neonic <- bind_rows(just_neonic, just_NWIS)
  
  return(nwis_neonic)
}
