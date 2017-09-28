library(readxl)
library(dplyr)
library(readr)

get_classes <- function(file_classes){
  
  classes <- data.frame(read_csv(file_classes))
  
  return(classes)
}

get_sites <- function(tracking,file.sites){
  dfSites <- as.data.frame(read_excel(file.sites))
  timeZone <- c("CST6CDT","CST6CDT","CST6CDT","EST5EDT","EST5EDT","EST5EDT","EST5EDT")
  States <- c("MN", "WI", "IN", "MI", "OH", "NY")
  names(timeZone) <- States
  
  names(dfSites) <- make.names(names(dfSites))
  dfSites$State <- substr(dfSites$Site.name,nchar(dfSites$Site.name)-1,nchar(dfSites$Site.name))
  dfSites$timeZone <- timeZone[dfSites$State]
  dfSites$shortName <- c("Bad","Manitowoc","IHC","St. Joe","Grand","Saginaw","Rouge","Maumee","Cuyahoga","Genesee")
  
  dfSites[which(dfSites$USGS.station.number == "04157005"),"USGS.station.number"] <- "04157000"
  
  from_NWIS <- readNWISsite(zeroPad(unique(tracking$SiteID),padTo = 8))
  
  dfSites <- left_join(from_NWIS, dfSites, by=c("site_no"="USGS.station.number"))
  
  return(dfSites)
}

get_neonic_data <- function(file.neonic,
                            file.MDL,
                            dfSites,
                            dfTracking){

  dfNeonic <- read.csv(file.neonic,stringsAsFactors = FALSE)
  
  dfMDL <- read.csv(file.MDL,stringsAsFactors = FALSE)
  
  #remove rows without sample information
  dfNeonic <- dfNeonic[grep("WS",dfNeonic$Sample,ignore.case = TRUE),]
  
  #Reconcile differences in site names to make them unique
  dfNeonic$Site <- sub(" at "," @ ",dfNeonic$Site)
  dfNeonic$Site[grep("Saginaw",dfNeonic$Site)] <- "Saginaw R @ Saginaw, MI" 
  dfNeonic$Site[grep("Cuyahoga",dfNeonic$Site)] <- "Cuyahoga R @ Saginaw, MI" 

  neonics <- names(dfNeonic)[11:16]
  
  # set NA values to MDL and create remark columns with "<" and estimated values
  for(i in 1:length(neonics)){
    MDL <- dfMDL[which(dfMDL$neonic==neonics[i]),"MDL"]
    naRows <- which(is.na(dfNeonic[,neonics[i]]))
    dfNeonic[naRows,neonics[i]] <- MDL
    dfNeonic$remarkCol <- NA
    dfNeonic[naRows,"remarkCol"] <- "<"
    estRows <- which(dfNeonic[,neonics[i]] < 0)
    dfNeonic[estRows,neonics[i]] <- -dfNeonic[estRows,neonics[i]]
    dfNeonic[estRows,"remarkCol"] <- "estimated"                 
    names(dfNeonic)[dim(dfNeonic)[2]] <- paste0("R_",neonics[i])
  }
  
  ## Merge tracking data with neonic results
  dfTracking <- as.data.frame(dfTracking)
  dfNeonic$Sample <- toupper(dfNeonic$Sample)
  dfTracking$Neonics <- toupper(dfTracking$Neonics)
  
  df <- dfNeonic %>%
    mutate(Date = as.Date(Date)) %>%
    left_join(select(dfTracking, -Site),
                  by = c("Sample"="Neonics")) %>%
    select(-Date.x, -Time.y, -Time.x, -Date.y) 

  return(df)
}