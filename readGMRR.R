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

gmrr <- bind_rows(gmrr1,gmrr2)
rm(gmrr1,gmrr2)

# analysis per cell
neededcell <- "D39773"
gMrrCell <- gmrr %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter (Cell == neededcell) %>%
  select(1:7,contains("RXLEV")) %>%
  gather("Measure","Value",8:135) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  separate(Measure,into=c("kpi","dir"),sep=5) %>%
  mutate(Bin=Bin-110)

ggplot(gMrrCell,aes(x=Bin,y=Value,fill=dir)) + geom_col(position = "dodge")

gMrrSmall <- gmrr %>%
  #filter(Band=="GSM900") %>%
  filter (NoReportsPassedFilter>1000) %>%
  select(1,5,7,338:369) %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  group_by(Cell,BSC) %>%
  summarize_all(funs(sum(.))) %>%
  ungroup() %>%
  unite(Ta,5:36,sep=",") %>%
  mutate(TaVec=strsplit(Ta,",")) %>%
  select(-Ta) %>%
  rowwise() %>%
  mutate(TA_50Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.5)[[1]]) %>%
  mutate(TA_75Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.75)[[1]]) %>%
  mutate(TA_90Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.9)[[1]])


MrrTa <- gMrrSmall %>%
  select (-NoReportsPassedFilter,-TaVec)

write_csv(MrrTa,"allDelta_Ta90Perc.csv")
