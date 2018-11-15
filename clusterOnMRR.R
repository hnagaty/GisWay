# I should run the file "readGenericMRR.R" before this
# the above script reads the names of complete MRR export files in a given folder

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(fastcluster)
library(tibble)
library(purrr)

myPath="d:/data/mrr/2018Sep25/"
mrrFiles=paste0(mrrConf.df$File,".msmt")
exportPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"

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
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter(NoReportsPassedFilter>1000) %>%
  select(-(NoFERULPassedFilter:NoFERULDLUnfiltered)) %>%
  filter(complete.cases(.))

write_csv(gmrr,paste0(exportPath,"GsmMRR.csv"))
gmrr <- read_csv(paste0(exportPath,"GsmMRR.csv"))
nCols <- ncol(gmrr)

gmrrRatio <- gmrr %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter)) %>%
  select(-contains("TrafficLevelM"))
gmrrRatio <- add_column(gmrrRatio, cluster=0, .after = 2)

data <- gmrrRatio %>%
  select(`TrafficLevelCS(E)`:`PATHLOSSDIFF(25,25)`)

dataScaled <- scale(data)
#mrr.dist <- dist(dataScaled, method = 'euclidean') #not needed as I use hclust,vector
#mrr.hclust.fit <- hclust(mrr.dist, method = 'complete') #causes memory overflow
mrrHclust.fit <- hclust.vector(dataScaled,method="ward",metric="euclidean")
saveRDS(mrrHclust.fit,file=paste0(exportPath,"HClustFit_Nov18.rds"))
mrrHclust.fit <- readRDS(file=paste0(exportPath,"HClustFit_Nov18.rds"))

# Dendogram plot
dendMrr <- as.dendrogram(mrrHclust.fit)
dendColored <- color_branches(dendMrr, k = 2)
plot(dendColored,leaflab="none")

mrr.clusters <- cutree(mrrHclust.fit, k = 3)
table(mrr.clusters)
gmrrRatio$cluster <- mrr.clusters
table(gmrrRatio$cluster,gmrrRatio$Band) # ==> cluster 1 is mostly 1800, cluster is mostly 900

# Then, plots for each class
neededCell <- '55062'
gmrrTidy <-  gmrrRatio %>% # should make group_by cluster first step, to speed it up
  #filter(Cell==neededCell | Cell==paste0("C",neededCell)) %>%
  gather("Measure","Value",`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`) %>%
  select(BSC:Band,Measure,Value) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  mutate(dir=substring(Measure,regexpr("UL|DL",Measure))) %>%
  separate(Measure,into=c("kpi"),sep="UL|DL") %>%
  group_by(cluster,kpi,dir,Bin) %>%
  dplyr::summarize(Value=mean(Value)) %>%
  rowwise() %>%
  mutate(dir = pasteDir(dir)) %>%
  ungroup() %>%
  mutate(Bin=ifelse(kpi=="RXLEV",Bin-110,Bin))

#saveRDS(gmrrTidy,file=paste0(exportPath,"gMrrTidy.rds"))
#gmrrTidy <- readRDS(file=paste0(exportPath,"gMrrTidy.rds"))

kpiListBi1 <- c("PATHLOSS","RXLEV")
kpiListBi2 <- c("RXQUAL")
kpiListUni1 <- c("PATHLOSSDIFF")
kpiListUni2 <- c("BSPOWER","MSPOWER","TAVAL")

# Define generic plotting functions
plotBiDirLine <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=dir)) +
    geom_line(position = "dodge",size=1) +
    facet_grid(cluster~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",color="Direction")
}

plotBiDirBar <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,fill=dir)) +
    geom_col(position = "dodge",col="grey") +
    facet_grid(cluster~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",fill="Direction")
}

plotUniDirLine <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=cluster)) +
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
    facet_grid(cluster~.) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",fill="Direction")
}

map(kpiListBi1,plotBiDirLine)
map(kpiListBi2,plotBiDirBar)
map(kpiListUni1,plotUniDirLine)
map(kpiListUni2,plotUniDirBar,8)

