#Add flow to Neonic data

library(dplyr)
library(dataRetrieval)
library(readxl)
library(USGSHydroTools)

#Define data directories:
raw.path <- "raw_data"
cached.path <- "cached_data"
cached.read <- "1_munge"
cached.save <- "2_process"

file.sites <- "GLRISiteCharacteristicsNeonicsForR.xlsx"
file.data <- "NeonicConcLandChar.rds"

dfSites <- as.data.frame(read_excel(file.path(raw.path,file.sites)))
dfNeonic <- readRDS(file.path(cached.path,cached.read,file.data))
names(dfSites) <- make.names(names(dfSites))

#Add leading zero to site IDs in neonic df.
dfNeonic$USGS.Site.ID <- ifelse(!is.na(dfNeonic$USGS.Site.ID), paste0("0",dfNeonic$USGS.Site.ID),NA)

sites <- dfSites$USGS.station.number
for(i in 1:length(sites)){
 site <- sites[i]
  i <- 1
  dfQ <- readNWISuv(dfSites[i,"USGS.station.number",],"00060", "2015-10-01","2016-09-30")
  dfDMQ <- readNWISdv(dfSites[i,"USGS.station.number",],"00060", "2015-10-01","2016-09-30")
  unique(dfDMQ$X_00060_00003_cd)
  

varnames <- "X_00060_00000"
test <- TSstats(df=dfQ,date="dateTime", varnames=varnames,dates = dfNeonic, starttime = "pdate",
        times=6,units="hrs",stats.return = c("mean","max","min","nearprev","difference"),
        subdatesvar = "USGS.Site.ID",subdatesvalue = site,out.varname = "Q_")
