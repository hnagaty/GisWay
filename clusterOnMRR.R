# I should run the file "readGenericMRR.R" before this
# the above script reads the names of complete MRR export files in a given folder

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(fastcluster)

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
gmrr <- read_csv(paste0(outDir,"GsmMRR.csv"))

gmrr <- gmrr %>%
  separate(CellName,into=c("BSC","Cell"),sep="/") %>%
  filter(NoReportsPassedFilter>1000) %>%
  select(-(NoFERULPassedFilter:NoFERULDLUnfiltered)) %>%
  filter(complete.cases(.))
nCols <- ncol(gmrr)
gmrrRatio <- gmrr %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter)) %>%
  select(-contains("TrafficLevelM"))
data <- gmrrRatio %>%
  select(`TrafficLevelCS(E)`:`PATHLOSSDIFF(25,25)`)

dataScaled <- scale(data)
#mrr.dist <- dist(dataScaled, method = 'euclidean') #not needed as I use hclust,vector
#mrr.hclust.fit <- hclust(mrr.dist, method = 'complete') #causes memory overflow

mrr.hclust.fit <- hclust.vector(dataScaled,method="ward",metric="euclidean")

mrr.clusters <- cutree(mrr.hclust.fit, k = 4)
table(mrr.clusters)

plot(mrr.hclust.fit)
