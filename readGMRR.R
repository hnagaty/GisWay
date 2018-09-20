library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)

pasteDir <- function(c) {
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}
  

myPath="D:/Optimisation/~InProgress/201806_GisFramework/ossData/"
mrrFile=c("Delta01_20180801.msmt","Delta02_20180801.msmt")

gmrr <- data.frame()
for (m in mrrFile) {
  gmrrtmp <- read_tsv(paste0(myPath,m))
  gmrr <- bind_rows(gmrr,gmrrtmp)
}
rm(gmrrtmp,m,mrrFile,myPath)

#Remove the FER as it's alawys 0
gmrr <- gmrr %>%
  select(-starts_with("FER"))

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

nCols <- ncol(gmrr)
gmrrFilter <- gmrr %>%
  filter(NoReportsPassedFilter>1000)

gmmrRatio <- gmrrFilter %>%
  mutate_at(16:nCols,funs(. / NoReportsPassedFilter))

gmrrSmall <- gmrrFilter %>%
  select(CellName:NoFERULDLUnfiltered) %>%
  separate(CellName,into=c("BSC","Cell"),sep="/")  

gmrrTidy <-  gmrr %>%
  filter(NoReportsPassedFilter>1000) %>%
  mutate_at(16:nCols,funs(. / NoReportsPassedFilter)) %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter(Cell==neededcell) %>%
  gather("Measure","Value",17:nCols) %>%
  select(BSC:Band,Measure,Value) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  mutate(dir=substring(Measure,regexpr("UL|DL",Measure))) %>%
  separate(Measure,into=c("kpi"),sep="UL|DL") %>%
  rowwise() %>%
  mutate(dir = pasteDir(dir)) %>%
  ungroup() %>%
  select(BSC:kpi,dir,Bin,Value)

test <-gmrrTidy %>%
  group_by(BSC,Cell,ChannelGroup,SubCellType,Band,kpi,dir) %>%
  select(-Bin) %>%
  dplyr::summarize(s=sum(Value,na.rm=TRUE))

test2 <- gmrrTidy %>%
  filter(kpi=="FER") %>%
    select(Value)


gmmrTidy <- gmrrTidy %>%
  mutate(TA_10Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.1)[[1]]) %>%
  mutate(TA_50Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.5)[[1]]) %>%
  mutate(TA_90Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.9)[[1]])

distinct(gmrrTidy,kpi)
distinct(gmrrTidy,Value)
a <- distinct(gmrrTidy,kpi,Bin)
write_csv(a,"ranges.csv")

gmrrTidy %>% filter(kpi=="RXQUAL") %>%
  ggplot(aes(x=Bin,y=Value,fill=dir)) + geom_col(position = "dodge") +
  facet_grid(Band~.)

gmrrTidy %>% filter(kpi=="PATHLOSSDIFF") %>%
  ggplot(aes(x=Bin,y=Value,col=dir)) +
  geom_jitter(alpha=0.4) +
  facet_grid(Band~dir)


gmrrTidy %>%
  filter(kpi=="RXQUAL",dir=="DL",Bin==61,Band=="GMS900",Value>100)


# I used this before, but now I want to make the one in the above lines
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
