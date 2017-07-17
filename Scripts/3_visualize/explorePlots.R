#Explore relations between neonic concentrations and other variables

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


plot(dfNeonic$Imidacloprid~jitter(dfNeonic$Ag..crops,factor = 20),log="y")
plot(dfNeonic$Clothianidin~jitter(dfNeonic$Ag..crops,factor = 20),log="y")
plot(dfNeonic$Thiamethoxam~jitter(dfNeonic$Ag..crops,factor = 20),log="y")
plot(dfNeonic$Thiamethoxam~dfNeonic$pdate,log="")

plot(dfNeonic$Thiamethoxam ~ dfNeonic$pdate)

library(smwrQW)
sinDate <-fourier(as.POSIXct(dfNeonic$Date,format="%m/%d/%Y"))[,1]
cosDate <-fourier(as.POSIXct(dfNeonic$Date,format="%m/%d/%Y"))[,2]

summary(lm(log10(dfNeonic$Clothianidin) ~dfNeonic$Ag..crops+sinDate+cosDate))
summary(lm(log10(dfNeonic$Thiamethoxam) ~dfNeonic$Ag..crops+sinDate+cosDate))


dfNeonic$pdate <- paste(dfNeonic$Date, dfNeonic$Time)
as.POSIXct(dfNeonic$pdate,format='%m/%d/%Y %H:%M',tz=timeZone)

