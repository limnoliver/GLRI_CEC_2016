library(readxl)
#library(data.table)
library(dplyr)
#library(dataRetrieval)
#library(lubridate)

#Raw data folder:
raw.path <- "raw_data"
cached.path <- "cached_data"
cached.read <- "0_download"
cached.save <- "1_munge"

timeZone <- c("CST6CDT","CST6CDT","CST6CDT","EST5EDT","EST5EDT","EST5EDT","EST5EDT")
States <- c("MN", "WI", "IN", "MI", "OH", "NY")
names(timeZone) <- States

#openState <- function(raw.path, cached.path, cached.save){
file.neonic <- "GLRI Neonic to NWIS 011317.csv"
file.sites <- "GLRISiteCharacteristicsNeonicsForR.xlsx"
file.tracking <- "tracking.rds"
file.MDL <- "neonicMDLs.csv"

dfNeonic <- read.csv(file.path(raw.path,file.neonic),stringsAsFactors = FALSE)
dfSites <- as.data.frame(read_excel(file.path(raw.path,file.sites)))
dfTracking <- readRDS(file.path(cached.path,cached.read,file.tracking))
dfMDL <- read.csv(file.path(raw.path,file.MDL),stringsAsFactors = FALSE)

#Extract state and determine time zones
names(dfSites) <- make.names(names(dfSites))
dfSites$State <- substr(dfSites$Site.name,nchar(dfSites$Site.name)-1,nchar(dfSites$Site.name))
dfSites$timeZone <- timeZone[dfSites$State]

#remove rows without sample information
dfNeonic <- dfNeonic[grep("WS",dfNeonic$Sample,ignore.case = TRUE),]

#Reconcile differences in site names to make them unique
dfNeonic$Site <- sub(" at "," @ ",dfNeonic$Site)
dfNeonic$Site[grep("Saginaw",dfNeonic$Site)] <- "Saginaw R @ Saginaw, MI" 
dfNeonic$Site[grep("Cuyahoga",dfNeonic$Site)] <- "Cuyahoga R @ Saginaw, MI" 
unique(dfNeonic$Site)

dfSites[which(dfSites$USGS.station.number == "04157005"),"USGS.station.number"] <- "04157000"

# #Make Neonic site names consistent with sites file
# siteNames <- character()
# siteIDs <- character()
# for(i in 1:dim(dfNeonic)[1]){
#   sitePre <- substr(dfNeonic[i,"Site"],1,5)
#   siteRow <- which(substr(dfSites$Site.name,1,5) == sitePre)
#   siteNames <- c(siteNames,dfSites[siteRow,"Site.name"])
#   siteIDs <- c(siteIDs,dfSites[siteRow,"USGS.station.number"])
# }

# dfNeonic$siteNames <- siteNames
# dfNeonic$siteIDs <- siteIDs


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
  


c(unique(siteNames))

## Merge tracking data with neonic results
dfNeonic$Sample <- toupper(dfNeonic$Sample)
dfTracking$Neonics <- toupper(dfTracking$Neonics)
df <- merge(dfNeonic,dfTracking,by.x = "Sample", by.y = "Neonics",all=TRUE)

#Merge watershed characteristics with Neonic data
as.numeric(dfSites$USGS.station.number) %in% dfNeonic$USGS.Site.ID
dfSiteInfo <- as.data.frame(dfNeonic$Sample)
for(i in 4:dim(dfSites)[2]){
  variable <- names(dfSites)[i]
  values <- dfSites[,variable]
  names(values) <- paste0("ID",dfSites$USGS.station.number)
  siteIDs <- ifelse(!is.na(dfNeonic$USGS.Site.ID), paste0("ID0",dfNeonic$USGS.Site.ID),NA)
  dfSiteInfo <- cbind(dfSiteInfo,values[siteIDs])
}
names(dfSiteInfo) <- c("Sample",names(dfSites)[-(1:3)])
dfNeonic <- merge(dfNeonic,dfSiteInfo,by="Sample",all=TRUE)
  

plot(dfNeonic$Imidacloprid~dfNeonic$Ag..crops)
plot(dfNeonic$Clothianidin~dfNeonic$Ag..crops,log="y")

summary(lm(dfNeonic$Imidacloprid~dfNeonic$Ag..crops+dfNeonic$Urban))


dfNeonic$pdate <- paste(dfNeonic$Date, dfNeonic$Time)
as.POSIXct(dfNeonic$pdate,format='%m/%d/%Y %H:%M',tz=timeZone)

