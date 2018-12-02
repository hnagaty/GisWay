# Reads MRR files
# Summarizes them
# Do some plots on idividual cells or groups


# Load libraries ----------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(purrr)


# Variables & functions ---------------------------------------------------
pasteDir <- function(c) {
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}

myPath <- "d:/data/mrr/2018Sep25/"
myPath <- "~/data/gmrr/"
myPath <- dataDir
mrrFiles <- paste0(mrrConf.df$File,".msmt")

outDir <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
#outDir <- "~/Dropbox/Voda/GISWay/export/"


# Read in the GSM MRR files -----------------------------------------------
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


# Save &  Load ------------------------------------------------------------
write_csv(gmrr,paste0(outDir,"GsmMRR.csv"))
gmrr <- read_csv(paste0(outDir,"GsmMRR.csv"))


# Further Processing ------------------------------------------------------
gmrrFiltered <- gmrr %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter(NoReportsPassedFilter>1000) %>%
  select(-(NoFERULPassedFilter:NoFERULDLUnfiltered)) %>%
  filter(complete.cases(.))
nCols <- ncol(gmrrFiltered)
gmrrRatio <- gmrrFiltered %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter))
#rm(gmrrFiltered,gmrr)

# Tidy for single cell plotting -------------------------------------------
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
  select(BSC:kpi,dir,Bin,Value) %>%
  mutate(Bin=ifelse(kpi=="RXLEV",Bin-110,Bin))


# KPI List ----------------------------------------------------------------
kpiListBi1 <- c("PATHLOSS","RXLEV")
kpiListBi2 <- c("RXQUAL")
kpiListUni1 <- c("PATHLOSSDIFF")
kpiListUni2 <- c("BSPOWER","MSPOWER","TAVAL")


# Plotting Functions ------------------------------------------------------
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
    mutate(Bin=as.factor(Bin)) %>%
    ggplot(aes(x=Bin,y=Value)) +
    geom_col(fill="orange") +
    facet_grid(Band~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",fill="Direction")
}


# Generate the Plots ------------------------------------------------------

map(kpiListBi1,plotBiDirLine)
map(kpiListBi2,plotBiDirBar)
map(kpiListUni1,plotUniDirLine)
map(kpiListUni2,plotUniDirBar,8)


