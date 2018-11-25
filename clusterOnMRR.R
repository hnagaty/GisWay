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
library(doParallel)
library(dendextend)
library(tictoc)


registerDoParallel(cores=2) #multi-core like functionality


pasteDir <- function(c) {
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}

dataDir<- "d:/data/mrr/2018Sep25/"
dataDir <- "~/data/gmrr/"
exportPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
exportPath <- "~/Dropbox/Voda/GISWay/export/"

mrrFiles=paste0(mrrConf.df$File,".msmt")
gmrr <- data.frame()
for (m in mrrFiles) {
  gmrrtmp <- read_tsv(paste0(dataDir,m))
  gmrr <- bind_rows(gmrr,gmrrtmp)
}
rm(gmrrtmp,m,mrrFiles)

#Remove the FER as it's alawys 0
gmrr <- gmrr %>%
  select(-starts_with("FER"))

#Combine all instances of a single cell
gmrr <- gmrr %>%
  group_by(CellName,ChannelGroup,SubCellType,Band) %>%
  summarise_all(sum) %>%
  ungroup() %>%
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

#rm(gmrr)

neededCell='D71194'

#filter(Cell==neededCell | Cell==paste0("C",neededCell)) %>%

qtls <- 0.5
qtlsnames <- "Median"
qtls <- seq(0.15,0.9,0.15)
qtlsnames <- paste0("Q",qtls*100)

gmrr2 <-top_n(gmrr,100)

tic("Tidy")
gmrrFeatures <-  gmrr2 %>%
  select(Cell,ChannelGroup,Band,`TrafficLevelCS(E)`,`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`) %>%
  rename(Traffic=`TrafficLevelCS(E)`) %>%  
  gather("Measure","Value",`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  group_by(Cell,ChannelGroup,Band,Measure) %>%
  summarise(traffic=max(Traffic),val=paste0(Bin,collapse=","),cnt=paste0(Value,collapse=",")) %>%
  rowwise() %>%
  mutate(q=paste0(wtd.quantile(as.numeric(strsplit(val,",",fixed=TRUE)[[1]]),
                                weights=as.numeric(strsplit(cnt,",",fixed=TRUE)[[1]]),
                                probs=qtls),
         collapse=",")) %>%
  select(-val,-cnt) %>%
  #separate(q,sep=",",into=qtlsnames)
  spread(key=Measure,value = q)
toc()

gmrrReshaled <- reshape(gmrrFeatures,idvar=c("Cell","ChannelGroup","Band","traffic"),timevar="Measure",direction = "wide")
saveRDS(gmrrTidy,file=paste0(exportPath,"gMrrTidy.rds"))
gmrrTidy <- readRDS(file=paste0(exportPath,"gMrrTidy.rds"))
  

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


data <- gmrrRatio %>%
  select(`TrafficLevelCS(E)`:`PATHLOSSDIFF(25,25)`)

dataScaled <- scale(data)
#mrr.dist <- dist(dataScaled, method = 'euclidean') #not needed as I use hclust,vector
#mrr.hclust.fit <- hclust(mrr.dist, method = 'complete') #causes memory overflow
mrrHclust.fit <- hclust.vector(dataScaled,method="ward",metric="euclidean")
saveRDS(mrrHclust.fit,file=paste0(exportPath,"HClustFit_Nov21.rds"))
mrrHclust.fit <- readRDS(file=paste0(exportPath,"HClustFit_Nov21.rds"))

# Dendogram plot
k=2
dendMrr <- as.dendrogram(mrrHclust.fit)
dendColored <- color_branches(dendMrr, k = k)
plot(dendColored,leaflab="none")

mrr.clusters <- as.factor(cutree(mrrHclust.fit, k = k))
table(mrr.clusters)
gmrrRatio$cluster <- mrr.clusters
table(gmrrRatio$cluster,gmrrRatio$Band) # ==> cluster 1 is mostly 1800, cluster is mostly 900

# Then, plots for each class
neededCell <- '55062'

rm(gmrrTidy)
tic("Tidy GMRR")
gmrrTidy <-  gmrrRatio %>%
  group_by(cluster) %>%
  summarise_at(vars(`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`),mean) %>%
  gather("Measure","Value",`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`) %>%
  separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
  select(-c,-d) %>%
  mutate(dir=substring(Measure,regexpr("UL|DL",Measure))) %>%
  separate(Measure,into=c("kpi"),sep="UL|DL",extra="drop") %>%
  rowwise() %>%
  mutate(dir = pasteDir(dir)) %>%
  ungroup() %>%
  mutate(Bin=ifelse(kpi=="RXLEV",Bin-110,Bin))
toc()


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
    geom_line(position = "identity",size=1) +
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

