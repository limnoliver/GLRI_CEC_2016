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
begin <- "2015-10-01"
end <- "2016-10-01"
timeSteps <- as.numeric((difftime(as.POSIXct(end),as.POSIXct(begin),units="days")+2)*96)
datesVector <- as.POSIXct(begin,tz="UTC") + 5*60*60 + seq(from = 1, to = timeSteps)*15*60
dfDates <- data.frame(dateTime=datesVector)
listQ <- list()

#Separate out each site, and add flow data for each
for(i in 1:length(sites)){
  site <- sites[i]
  dfQ <- readNWISuv(dfSites[i,"USGS.station.number",],"00060", begin,end)
  dfDMQ <- readNWISdv(dfSites[i,"USGS.station.number",],"00060", begin,end)
  dfQ1 <- right_join(dfQ,dfDates,by = "dateTime")
  dfQ1$Date <- as.Date(dfQ1$dateTime)
  dfQ2 <- left_join(dfQ1,dfDMQ,by="Date")
  dfQ2$Q <- ifelse(!is.na(dfQ2$X_00060_00000),dfQ2$X_00060_00000,dfQ2$X_00060_00003)
  if(site == "04157005") site <- "04157000"
  subdfNeonic <- subset(dfNeonic,USGS.Site.ID==site)
  subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                  times=6,units="hrs",stats.return = c("mean","max","min","difference"),
                  subdatesvar = "USGS.Site.ID",subdatesvalue = site,out.varname = "Q")
  subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                         times=0,units="hrs",stats.return = c("nearprev"),
                         subdatesvar = "USGS.Site.ID",subdatesvalue = site,out.varname = "Q")
  subdfNeonic <- TSstats(df=dfQ2,date="dateTime", varnames="Q",dates = subdfNeonic, starttime = "pdate",
                         times=1,units="hrs",stats.return = c("mean"),
                         subdatesvar = "USGS.Site.ID",subdatesvalue = site,out.varname = "Q")
  
  listQ[[i]] <- subdfNeonic
}

#Combine individual site data frames back into one with all flow data included
dfNeonic2 <- listQ[[1]]
for(i in 2:length(sites)){
  dfNeonic2 <- cbind(dfNeonic2,listQ[[i]])
}
                     

plot(test2$Q~test2$dateTime,lty=1,pch=".")
plot(test2$X_00060_00003~test2$dateTime,lty=1,pch=".")

dateValues <- as.data.frame(table(as.Date(dfQ$dateTime)))
unique(dfDMQ$X_00060_00003_cd)



varnames <- "X_00060_00000"
