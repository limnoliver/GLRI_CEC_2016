# Read data from Google Tracking:
source("scripts/0_download/getGoogleSampleTrackingDownload.R")
#' @examples 
#' tracking <- readRDS(file.path("cached_data","0_download","tracking.rds"))

# Open neonic data:
source("scripts/1_munge/openStateData.R")
#' @examples 
#' data.mi <- readRDS(file.path("cached_data","1_munge","dataMI.rds"))
#' data.wi <- readRDS(file.path("cached_data","1_munge","dataWI.rds"))
#' data.ny <- readRDS(file.path("cached_data","1_munge","dataNY.rds"))
# total rows 471, which includes blanks and QAs

# Merge GLRI/GLPF optic fluorescence/absorbance data:
# GLRIWWOptFlAbsVectorized.RData
source("scripts/2_process/mergeGLPF_GLRI_opt.R", encoding = 'UTF-8')
#' @examples 
#' dfabs <- readRDS(file.path("cached_data","2_process_merge_GLRI_GLPF","dfabs.rds"))
#' dffl <- readRDS(file.path("cached_data","2_process_merge_GLRI_GLPF","dffl.rds"))
# This is all optics data, QA and blanks

# Merge tracking and bacteria data
# Files that make me nervous because I have no idea where they came from:
# GLRI01-04-16_mergedBact.RData
# GLRIWWMar162016summary.RData
source("scripts/2_process/mergeTrackingBacteriaData.R", encoding = 'UTF-8')
#' @examples 
#' trackingBacteria <- readRDS(file.path("cached_data","3_process_merge_bacteria","trackingBacteria.rds"))

# Merge field and opt data
source("scripts/2_process/mergeFieldOpt.R", encoding = 'UTF-8')
#' @examples
#'mergedData <- readRDS(file.path("cached_data","4_process_merge_field_opt","mergedData.rds"))

# Separate sanitary, ww, qa, blanks, etc
source("scripts/2_process/filterData.R", encoding = 'UTF-8')
#' @examples
#' mergedData <- readRDS(file.path("cached_data","5_process_filterData","rds","summary_noWW_noQA.rds"))
#' only.data.with.field.data <- mergedData[rowSums(is.na(mergedData[,c("WT","DO","Turb","SC","pH")])) != 5,]

source("scripts/2_process/categorizeBacteria.R", encoding = 'UTF-8')
#' @export
#' mergedData <- readRDS(file.path("cached_data","6_process_cagegorize_Bacteria","rds","summary_noWW_noQA.rds"))

source("scripts/2_process/summaryOpticalVariables.R", encoding = 'UTF-8')
#' @examples
#' summaryDF <- readRDS(file.path("cached_data","7_process_summarize_optics","rds",paste0("summary","_noQA",".rds")))


source("scripts/2_process/newCategories.R", encoding = 'UTF-8')
#' @examples
#' summaryDF <- readRDS(file.path("cached_data","8_process_new_categories","rds",paste0("summary","_noQA",".rds")))

#########################################################
## ADD EventFreqAndDatesDetermination.R into workflow  
## The file name was changed after running that script 
## and an event grouping was added by hand
#########################################################

source("scripts/2_process/eventGrouping.R", encoding = 'UTF-8')
#' @examples
#' summaryDF <- readRDS(file.path("cached_data","8_process_new_categories","rds",paste0("summary","_noQA",".rds")))


summaryDF <- readRDS(file.path("cached_data","8_process_new_categories","rds",paste0("summary","_noQA",".rds")))
dfabs <- readRDS(file.path("cached_data","5_process_filterData","rds","dfabs_noWW_noQA.rds"))
dffl <- readRDS(file.path("cached_data","5_process_filterData","rds","dffl_noWW_noQA.rds"))
EEMS <- readRDS(file.path("cached_data","7_process_summarize_optics","rds","EEMs3D_noWW_noQA.rds"))

na.info <- function(df, key = "CAGRnumber", first.col = "OB1"){
  key.index <- which(names(df) == key)
  opt.df <- df[,c(key.index,which(names(df) == first.col):ncol(df))]
  df.noNA <- na.omit(opt.df)
  df.NA <- opt.df[!(opt.df[[key]] %in% df.noNA[[key]]),]
  na.cols.full <- names(opt.df)[!(names(opt.df) %in% names(df.noNA))]
  na.cols.partial <- colnames(df.NA)[ apply(df.NA, 2, anyNA) ]
  na.rows <- df.NA[[key]]
  
  inf.cols <- names(opt.df)[unlist(do.call(data.frame,lapply(opt.df,
                                                             function(x) any(is.infinite(x)))))]
  inf.rows <- which(is.infinite(rowSums(opt.df[-1])))
  
  nan.cols <- names(opt.df)[unlist(do.call(data.frame,lapply(opt.df,
                                                             function(x) any(is.nan(x)))))]
  nan.rows <- which(is.nan(rowSums(opt.df[-1])))
  
  return(list(na.cols.full = na.cols.full,
              na.cols.partial = na.cols.partial,
              na.rows = na.rows,
              inf.cols = inf.cols,
              inf.rows = inf.rows,
              nan.cols = nan.cols,
              nan.rows = nan.rows))
}


na.info.list <- na.info(summaryDF)

rmRows <- unique(c(which(summaryDF$CAGRnumber %in% na.info.list$na.rows),
                   na.info.list$nan.rows,
                   na.info.list$inf.rows))
summaryDF <- summaryDF[-rmRows,]

df.NA <- summaryDF %>% 
  select(c(1,match(na.info.list$na.cols.partial,names(.)))) %>%
  filter(CAGRnumber %in% na.info.list$na.rows)

df.noNA <- summaryDF %>% 
  select(-match(na.info.list$na.cols.partial,names(.))) %>%
  select(-match(na.info.list$inf.cols,names(.))) %>%
  filter(!(CAGRnumber %in% na.info.list$na.rows))

###################################
# Graphing functions

source("scripts/3_visualize/summarizeOptics.R", encoding = 'UTF-8')
#' @examples
#' shell.exec(file.path("cached_figures",paste0("EEMs","_noWW_noQA",".pdf")))
#' shell.exec(file.path("cached_figures",paste0("Abs","_noWW_noQA",".pdf")))
#' 

source("scripts/4_model/regressionTrees.R", encoding = 'UTF-8')
