library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)

#myPath="D:/Optimisation/~InProgress/201806_GisFramework/ossData/"
#mrrFile="Delta01_20180801.msmt"
#gmrr <- read_tsv(paste0(c(myPath,mrrFile)))

gmrr1 <- read_tsv("D:/Optimisation/~InProgress/201806_GisFramework/ossData/Delta01_20180801.msmt")
gmrr2 <- read_tsv("D:/Optimisation/~InProgress/201806_GisFramework/ossData/Delta02_20180801.msmt")

sd <- gmrr2 %>% separate(CellName,into=c("BSC","Cell"),sep="/") 
unique(sd$BSC)

gMrrSmall <- gmrr %>%
  filter(Band=="GSM900", NoReportsPassedFilter>1000) %>%
  select(1,5,7,338:369) %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  group_by(Cell,BSC) %>%
  summarize_all(funs(sum(.))) %>%
  ungroup() %>%
  unite(Ta,5:36,sep=",",) %>%
  mutate(TaVec=strsplit(Ta,",")) %>%
  select(-Ta) %>%
  rowwise() %>%
  mutate(Ta90Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.9)[[1]])

MrrTa <- gMrrSmall %>%
  select (Cell,Ta90Perc)

write_csv(MrrTa,"ta.csv")
