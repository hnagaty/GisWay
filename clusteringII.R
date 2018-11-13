# clustering trials
# based on ho
# later make this with mrrs
# and then check for multi-variate outliers

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

ossDir <- "~/Dropbox/Voda/GISWay/ossData/"
hoFile <- "GsmHoGrouped.csv"

hoRaw <- read_csv(paste0(ossDir,hoFile),
                  col_types = cols(
                    CellName = col_character(),
                    NbrCellName = col_character(),
                    AdjType = col_character(),
                    Daytime = col_character(),
                    DayType = col_character(),
                    HoCnt = col_double(),
                    HoSucc = col_double(),
                    CellTraffic = col_double()
                  ))
hoSpread <- hoRaw %>%
  select(-HoSucc,-CellTraffic) %>%
  filter(HoCnt!=0) %>%
  unite(Period,DayType,Daytime,sep="_") %>%
  unite(Relation,CellName,NbrCellName,sep="-") %>%
  spread(Period,HoCnt,fill=0)
