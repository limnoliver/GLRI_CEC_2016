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

  df <- merge(dfNeonic,dfTracking,by.x = "Sample", by.y = "Neonics",all=TRUE)
  df <- df[!is.na(df$Clothianidin),]
  
  #Populate Neonic file with land use characteristics and other site info
  dfNeonicSiteIDs <- ifelse(!is.na(dfNeonic$USGS.Site.ID), paste0("0",dfNeonic$USGS.Site.ID),NA)
  siteRows <- match(dfNeonicSiteIDs,dfSites$USGS.station.number)
  dfSiteInfo <- dfSites[siteRows,]
  dfNeonic <- cbind(dfNeonic,dfSiteInfo)
  
  timeZone <- c("CST6CDT","CST6CDT","CST6CDT","EST5EDT","EST5EDT","EST5EDT","EST5EDT")
  States <- c("MN", "WI", "IN", "MI", "OH", "NY")
  names(timeZone) <- States
  
  #Convert dates and times to POSIXct in GMT
  dfNeonic$timeZone <- ifelse(is.na(dfNeonic$timeZone),"EST5EDT",dfNeonic$timeZone)
  uniqueTimeZones <- unique(timeZone)
  listState <- list()
  for(i in 1:length(uniqueTimeZones)){
    dfState <- subset(dfNeonic,timeZone==uniqueTimeZones[i])
    dfState$pdate <- as.POSIXct(paste(dfState$Date, dfState$Time),format='%m/%d/%Y %H:%M',tz=uniqueTimeZones[i])
    dfState$pdate <- as.POSIXct(format(as.POSIXct(dfState$pdate),tz="UTC",usetz=TRUE),tz="UTC")
    listState[[i]] <- dfState
  }
  
  dfNeonic <- rbind(listState[[1]],listState[[2]])
  
  return(dfNeonic)
}