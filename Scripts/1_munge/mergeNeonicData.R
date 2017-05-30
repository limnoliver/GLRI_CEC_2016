library(readxl)
#library(data.table)
library(dplyr)
#library(dataRetrieval)
#library(lubridate)

#Raw data folder:
raw.path <- "raw_data"
cached.path <- "cached_data"
cached.save <- "1_munge"

timeZone <- c("CST6CDT","CST6CDT","CST6CDT","EST5EDT","EST5EDT","EST5EDT","EST5EDT")
States <- c("MN", "WI", "IN", "MI", "OH", "NY")
names(timeZone) <- States

#openState <- function(raw.path, cached.path, cached.save){
  file.neonic <- "GLRI Neonic Summary 090216.csv"
  file.sites <- "GLRISiteCharacteristicsNeonicsForR.xlsx"
  
  dfNeonic <- read.csv(file.path(raw.path,file.neonic),stringsAsFactors = FALSE)
  dfSites <- read_excel(file.path(raw.path,file.sites))
  
  #Extract state and determine time zones
  names(dfSites) <- make.names(names(dfSites))
  dfSites$State <- substr(dfSites$Site.name,nchar(dfSites$Site.name)-1,nchar(dfSites$Site.name))
  dfSites$timeZone <- timeZone[dfSites$State]
  
  #remove rows without sample information
  dfNeonic <- dfNeonic[grep("WS",dfNeonic$Sample),]
  
  #Reconcile differences in site names to make them unique
  dfNeonic$Site <- sub(" at "," @ ",dfNeonic$Site)
  dfNeonic$Site[grep("Saginaw",dfNeonic$Site)] <- "Saginaw R @ Saginaw, MI" 
  dfNeonic$Site[grep("Cuyahoga",dfNeonic$Site)] <- "Cuyahoga R @ Saginaw, MI" 
  unique(dfNeonic$Site)
  
  #Make Neonic site names consistent with sites file
  siteNames <- character()
  siteIDs <- character()
  for(i in 1:dim(dfNeonic)[1]){
    sitePre <- substr(dfNeonic[i,"Site"],1,5)
    siteRow <- which(substr(dfSites$Site.name,1,5) == sitePre)
    siteNames <- c(siteNames,dfSites[siteRow,"Site.name"])
    siteIDs <- c(siteIDs,dfSites[siteRow,"USGS.station.number"])
  }
  
  dfNeonic$siteNames <- siteNames
  dfNeonic$siteIDs <- siteIDs
  
  
  ## NEED TO MERGE THE SITE INFO WITH THE NEONIC DATA
  
  
  grep("Maumee",unique(c(dfNeonic$Site,dfSites$Site.name)))
  
  
  dfNeonic$pdate <- paste(dfNeonic$Date, dfNeonic$Time)
  as.POSIXct(dfNeonic$pdate,format='%m/%d/%Y %H:%M',tz=timeZone)
