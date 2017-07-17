#Add flow to Neonic data

library(dplyr)
library(dataRetrieval)

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
dfSites$'USGS station number'

i <- 1
test <- readNWISuv(dfSites[i,"USGS.station.number",],"00060", "2015-10-01","2016-09-30")

?readNWISuv
