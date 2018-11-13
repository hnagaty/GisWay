# Reads MRR files
# Summarizes them
# Do some plots on idividual cells or non groups
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(purrr)

pasteDir <- function(c) {
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}
  

myPath="d:/data/mrr/2018Sep25/"
mrrFiles=paste0(mrrConf.df$File,".msmt")
outDir <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"

gmrr <- data.frame()
for (m in mrrFiles) {
  gmrrtmp <- read_tsv(paste0(myPath,m))
  gmrr <- bind_rows(gmrr,gmrrtmp)
}
rm(gmrrtmp,m,mrrFiles,myPath)

#Remove the FER as it's alawys 0
gmrr <- gmrr %>%
  select(-starts_with("FER"))

#Combine all instances of a single cell
gmrr <- gmrr %>%
  group_by(CellName,ChannelGroup,SubCellType,Band) %>%
  summarise_all(sum) %>%
  ungroup()

write_csv(gmrr,paste0(outDir,"GsmMRR.csv"))

gmrrFiltered <- gmrr %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter(NoReportsPassedFilter>1000) %>%
  select(-(NoFERULPassedFilter:NoFERULDLUnfiltered)) %>%
  filter(complete.cases(.))
nCols <- ncol(gmrrFiltered)
gmrrRatio <- gmrrFiltered %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter))

# Make a tidy version (from cols to rows) for the MRR
# For use in plotting
# It's okay for the warning messages below
# Too slow if not filtered by cellname

neededCell <- '55062'
gmrrTidy <-  gmrrRatio %>%
  filter(Cell==neededCell | Cell==paste0("C",neededCell)) %>%
  gather("Measure","Value",`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`) %>%
  select(BSC:Band,Measure,Value) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  mutate(dir=substring(Measure,regexpr("UL|DL",Measure))) %>%
  separate(Measure,into=c("kpi"),sep="UL|DL") %>%
  rowwise() %>%
  mutate(dir = pasteDir(dir)) %>%
  ungroup() %>%
  select(BSC:kpi,dir,Bin,Value)

kpiListBi1 <- c("PATHLOSS","RXLEV")
kpiListBi2 <- c("RXQUAL")
kpiListUni1 <- c("PATHLOSSDIFF")
kpiListUni2 <- c("BSPOWER","MSPOWER","TAVAL")

# Define generic plotting functions
plotBiDirLine <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=dir)) +
    geom_line(position = "dodge",size=1) +
    facet_grid(Band~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",color="Direction")
}

plotBiDirBar <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,fill=dir)) +
    geom_col(position = "dodge",col="grey") +
    facet_grid(Band~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",fill="Direction")
}

plotUniDirLine <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=Band)) +
    geom_line(position = "dodge",size=1) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",color="Direction")
}

plotUniDirBar <- function(kpiV,taLimit) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    filter(!(kpi=="TAVAL" & Bin >= taLimit )) %>%
    ggplot(aes(x=Bin,y=Value)) +
    geom_col(fill="orange") +
    facet_grid(Band~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",fill="Direction")
}

map(kpiListBi1,plotBiDirLine)
map(kpiListBi2,plotBiDirBar)
map(kpiListUni1,plotUniDirLine)
map(kpiListUni2,plotUniDirBar,8)

# Ignore below code chunk for now
#========================================================================
gmmrTidy <- gmrrTidy %>%
  mutate(TA_10Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.1)[[1]]) %>%
  mutate(TA_50Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.5)[[1]]) %>%
  mutate(TA_90Perc=wtd.quantile(seq(0,31,1),weights=as.numeric(unlist(TaVec)),probs=0.9)[[1]])

distinct(gmrrTidy,kpi)
distinct(gmrrTidy,Value)
a <- distinct(gmrrTidy,kpi,Bin)
write_csv(a,"ranges.csv")

gmrrTidy %>%
  filter(kpi=="RXQUAL",dir=="DL",Bin==61,Band=="GMS900",Value>100)
#========================================================================



# Calculate the percentile values
# Needs to be verified. Compare it with output of TA line
gmrrPercentile <- gmrrTidy %>%
  #filter(Cell==neededCell) %>%  
  arrange (BSC,Cell,ChannelGroup,kpi,dir,Bin) %>%
  group_by(BSC,Cell,ChannelGroup,SubCellType,Band,kpi,dir) %>%
  mutate(pcnt=Value/sum(Value),cumPcnt=cumsum(pcnt)) %>%
  select(-Value,-pcnt) %>%
  arrange(BSC,Cell,ChannelGroup,SubCellType,Band,kpi,dir,desc(cumPcnt)) %>%
  filter(cumPcnt<0.9) %>%
  top_n(1,cumPcnt) %>% top_n(1,Bin) %>%
  select(-cumPcnt) %>%
  unite(KPI,kpi,dir) %>%
  spread(KPI,Bin)

write_csv(gmrrPercentile,"GsmMrr.csv")

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




# analysis per cell
neededcell <- "11181"
gMrrCell <- gmrr %>%
  #separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter (Cell == neededcell) %>%
  select(1:7,contains("RXLEV")) %>%
  gather("Measure","Value",8:135) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  separate(Measure,into=c("kpi","dir"),sep=5) %>%
  mutate(Bin=Bin-110)

ggplot(gMrrCell,aes(x=Bin,y=Value,fill=dir)) + geom_col(position = "dodge")
