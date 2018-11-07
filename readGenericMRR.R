library(readr)
library(dplyr)

dataDir <- "d:/data/mrr/2018Sep25/"
allMrrFiles <- list.files(dataDir,pattern=".conf")
allMrrFiles <- sub(".conf","",allMrrFiles)

mrrConf.df <- data.frame(File=as.character(),
                         Creator=as.character(),
                         StartTime=as.character(),
                         EndTime=as.character(),
                         noCells=as.numeric(),
                         Exists=as.logical())

for (mrrFile in allMrrFiles) {
  mrrConf <- read_lines(paste0(dataDir,mrrFile,".conf"))
  mrrStartTime <-  mrrConf[grepl("start time", tolower(mrrConf))]
  mrrStartTime <- substring(mrrStartTime,regexpr("\t",mrrStartTime)+1)
  mrrEndTime <-  mrrConf[grepl("stop time", tolower(mrrConf))]
  mrrEndTime <- substring(mrrEndTime,regexpr("\t",mrrEndTime)+1)
  mrrCreator <-  mrrConf[grepl("creator", tolower(mrrConf))]
  mrrCreator <- substring(mrrCreator,regexpr("\t",mrrCreator)+1)
  mrrNoCells <-  mrrConf[grepl("no of cells", tolower(mrrConf))]
  mrrNoCells <- as.numeric(substring(mrrNoCells,regexpr("\t",mrrNoCells)+1))
  mrrExists <- file.exists(paste0(dataDir,mrrFile,".msmt"))
  mrr.df <- data.frame(File=mrrFile,
                      Creator=mrrCreator,
                      StartTime=mrrStartTime,
                      EndTime=mrrEndTime,
                      noCells=mrrNoCells,
                      Exists=mrrExists)
  mrrConf.df <- rbind.data.frame(mrrConf.df,mrr.df)
}

rm(mrr.df)
rm(mrrFile,mrrCreator,mrrStartTime,mrrEndTime,mrrNoCells,mrrExists)
rm(mrrConf)

mrrConf.df$StartTime <- as.POSIXct(mrrConf.df$StartTime)
mrrConf.df$EndTime <- as.POSIXct(mrrConf.df$EndTime)

mrrConf.df <- mrrConf.df %>%
  filter(Exists=="TRUE") %>%
  arrange(desc(StartTime))

