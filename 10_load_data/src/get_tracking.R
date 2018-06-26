# Need to deal with times with no colon

library(httr)
library(googlesheets)
library(dplyr)
library(yaml)

# Make sure Bison Connect is your default, and log out for good measure any personal gmail

get_tracking_data <- function(get_tracking.config = "10_load_data/cfg/tracking_config.yaml"){
  
  config.args <- yaml.load_file(get_tracking.config)
  
  set_config(config(ssl_verifypeer = 0L))
  token <- gs_auth(cache = FALSE)
  my_sheets <- gs_ls()
  
  glpfTitle <- gs_title(config.args$fetch.args[["gs_title"]])
  # 
  timeZone <- c("CST6CDT","CST6CDT","CST6CDT","EST5EDT","EST5EDT","EST5EDT","EST5EDT")
  States <- c("MN", "WI", "IN", "MI", "OH", "NY")
  names(timeZone) <- States
  
  for (i in 1:length(States)){
    dfState <- gs_read(glpfTitle,ws=States[i],range="A4:P100")
    names(dfState) <- c("Site","SiteID","Date","Time", "MediumCode","HydrologicCondition",
                        "SampleTypeCode","NWISRecordNumber",
                        "Glyphosate","Neonics","ToxCast","Metabolomics","Transcriptomics","Passive",
                        "Comments","GlyphosateComments")
    
    dfState <- dfState[which(!is.na(dfState$Date)),]
    dfState$SampleTypeCode <- as.character(dfState$SampleTypeCode)
    dfState$NWISRecordNumber <- as.character(dfState$NWISRecordNumber)
    if(class(dfState$Time)[1] %in% c("character","integer")){
      dfState$Time <- as.integer(dfState$Time)
      dfState$Time <- substr(as.POSIXct(sprintf("%04.0f", dfState$Time), format='%H%M'), 12, 16)
        } else dfState$Time <- format(strptime(dfState$Time, format="%H:%M"), format = "%H:%M")
    dfState$Time  <- as.character(dfState$Time)
    dfState$pdate <- paste(dfState$Date,dfState$Time)
    dfState$pdate <- as.POSIXct(dfState$pdate,format='%m/%d/%Y %H:%M',tz=timeZone[i])
    dfState$pdate <- as.POSIXct(format(as.POSIXct(dfState$pdate),tz="GMT",usetz=TRUE),tz="GMT")
    dfState$Date <- as.Date(dfState$Date,format='%m/%d/%Y' )
    # 
    dfState$State <- States[i]
    dfState <- subset(dfState, select=c(Site:Time,pdate,MediumCode:State))
    dfState <- subset(dfState, select=c(State,Site:GlyphosateComments))
    
    if(i == 1) {
      df <- dfState
    } else {
      df <- bind_rows(df,dfState)
    }
  }

  df$NWISRecordNumber <- gsub(" ", "", df$NWISRecordNumber)
  df$NWISRecordNumber <- gsub("(DB04)", "", df$NWISRecordNumber, fixed=TRUE)
  df$NWISRecordNumber <- gsub("(DB01)", "", df$NWISRecordNumber, fixed=TRUE)
  df$NWISRecordNumber <- gsub("(DB01", "", df$NWISRecordNumber, fixed=TRUE)
  df$NWISRecordNumber <- gsub("q", "", df$NWISRecordNumber)
  df$NWISRecordNumber <- gsub("(DB4)", "", df$NWISRecordNumber, fixed=TRUE)
  df$NWISRecordNumber <- gsub("(DB3)", "", df$NWISRecordNumber, fixed=TRUE)
  df$NWISRecordNumber <- zeroPad(df$NWISRecordNumber, 8)
  
  tracking <- df
  
  # Making things right:
  tracking$SiteID[tracking$SiteID == "04157005"] <- "04157000"
  #tracking <- filter(tracking, SampleTypeCode == "9")
  #tracking$Date[tracking$NWISRecordNumber == "01601472"] <- as.Date("2016-04-07")
  tracking$Date[tracking$NWISRecordNumber == "01600121"] <- as.Date("2016-02-02")
  tracking$SiteID[tracking$SiteID == "4249000"] <- "04249000"
  tracking$Time[tracking$NWISRecordNumber == "01600388"] <- "10:20"
  tracking$Time[tracking$NWISRecordNumber == "01600462"] <- "09:10"
  tracking$Date[tracking$NWISRecordNumber == "01600384"] <- as.Date("2015-12-09")
  tracking$Date[tracking$NWISRecordNumber == "01600381"] <- as.Date("2015-12-08")
  tracking$Time[tracking$NWISRecordNumber == "01601976"] <-  "08:00"
  tracking$Date[tracking$NWISRecordNumber == "01600102"] <- as.Date("2015-11-12")
  tracking$Time[tracking$NWISRecordNumber == "01600461"] <- "12:10"
  tracking$Time[tracking$NWISRecordNumber == "016002001"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01600207"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01600236"] <- "12:00"
  tracking$Time[tracking$NWISRecordNumber == "01600484"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01600973"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01601028"] <- "12:00"
  #tracking$Date[tracking$NWISRecordNumber == "01601028"] <- as.Date("2016-05-10")
  tracking$Time[tracking$NWISRecordNumber == "01601556"] <- "11:00"
  tracking$Time[tracking$NWISRecordNumber == "01601641"] <- "10:00"
  tracking$Time[tracking$NWISRecordNumber == "01600618"] <- "13:00"
  tracking$Time[tracking$NWISRecordNumber == "01600700"] <- "10:00"
  
  tracking$pdate_new <- paste(tracking$Date,tracking$Time)
  
  df <- data.frame()
  for(st in States){
    track_st <- filter(tracking, State == st)
    track_st$pdate_1 <- as.POSIXct(track_st$pdate_new,format='%Y-%m-%d %H:%M',tz=timeZone[st])
    attr(track_st$pdate_1, "tzone") <- "UTC"
    df <- bind_rows(df, track_st)
  }
  
  df <- select(df, -pdate) %>%
    rename(pdate = pdate_1)
  
  #Getting the record numbers to line up in WQP:
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "01600425"] <- "01600815"		
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "016002001"] <- "01600123"		
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "01600236"] <- "01600262"		
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "01600484"] <- "01600487"		
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "01600700"] <- "01600698"		
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "01600769"] <- "01600768"		
  # tracking$NWISRecordNumber[tracking$NWISRecordNumber == "01601616"] <- "01601615"		
  # tracking$SiteID[tracking$SiteID == "04157005"] <- "04157000"
  
  return(df)
  
}

