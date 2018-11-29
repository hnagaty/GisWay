# I should run the file "readGenericMRR.R" before this
# the above script reads the names of complete MRR export files in a given folder

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(Hmisc)
library(fastcluster)
library(dendextend)
library(tictoc)
library(dummies)
library(sp)
library(rgdal)
library(spdplyr)
library(tmap)
library(RColorBrewer)
library(viridisLite)

# Parallel Processing -----------------------------------------------------
registerDoParallel(cores=2) #multi-core like functionality


# Functions & Indentifiers ------------------------------------------------
pasteDir <- function(c) {
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}

dataDir<- "d:/data/mrr/2018Sep25/"
#dataDir <- "~/data/gmrr/"
exportPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
#exportPath <- "~/Dropbox/Voda/GISWay/export/"
geoPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export"

#Plotting functions
kpiListBi1 <- c("PATHLOSS","RXLEV")
kpiListBi2 <- c("RXQUAL")
kpiListUni1 <- c("PATHLOSSDIFF")
kpiListUni2 <- c("BSPOWER","MSPOWER","TAVAL")

# Define generic plotting functions
plotBiDirLine <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=dir)) +
    geom_line(position = "dodge",size=1) +
    facet_wrap(~cluster,ncol=3) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",color="Direction")
}

plotBiDirBar <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,fill=dir)) +
    geom_col(position = "dodge",col="grey") +
    facet_wrap(~cluster,ncol=3) +
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

plotUniDirLineC <- function(kpiV) {
  gmrrTidy %>% filter(kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=cluster)) +
    facet_wrap(~cluster,ncol=3) +
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
    facet_wrap(~cluster,ncol=3) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle="Distribution",x="Value",y="Percentage",fill="Direction")
}



# Read MRR files ----------------------------------------------------------
# mrrconf.df should be pre-defined. It's defined in "readGenericMRR.R"
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

# Reduce Features --------------------------------------------------------------
gmrr <- read_csv(paste0(exportPath,"GsmMRR.csv"))
nCols <- ncol(gmrr)

gmrrRatio <- gmrr %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter)) %>%
  select(-contains("TrafficLevelM"))
#gmrrRatio <- add_column(gmrrRatio, cluster=0, .after = 3)

#rm(gmrr)

#neededCell='D71194'
#filter(Cell==neededCell | Cell==paste0("C",neededCell)) %>%

qtls <- 0.5
qtlsnames <- "Median"
qtls <- seq(0.15,0.9,0.15)
qtlsnames <- paste0("Q",qtls*100)

#gmrr2 <-top_n(gmrr,100)

tic("MRR Features")
gmrrFeatures <-  gmrr %>%
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
  separate(q,sep=",",into=qtlsnames,convert=TRUE) %>%
  gather(key=qtl,value=value,-(Cell:traffic)) %>%
  unite(ftr,Measure,qtl,sep="-") %>%
  spread(key=ftr,value = value) %>%
  ungroup()
toc()
gmrrFeatures$Band <- as.factor(gmrrFeatures$Band)

saveRDS(gmrrFeatures,file=paste0(exportPath,"gmrrFeatures.rds"))


# Clustering --------------------------------------------------------------
gmrrFeatures <- readRDS(file=paste0(exportPath,"gmrrFeatures.rds"))

data <- gmrrRatio %>%
  select(`TrafficLevelCS(E)`:`PATHLOSSDIFF(25,25)`)

removeFeatures <- vars(Band,traffic)
data <- gmrrFeatures %>%
  filter(Band=="GSM900") %>%
  select(-Cell,-ChannelGroup,-Band,-traffic)
#data <- select(gmrrFeatures,-!!!removeFeatures) # I want to work this out
#bandDummy <- dummy(gmrrFeatures$Band) #one hot encoding
#rownames(bandDummy) <- gmrrFeatures$Cell #Not actually needed
#data <- bind_cols(data,as.data.frame(bandDummy))
labels <- gmrrFeatures %>%
  filter(Band=="GSM900") %>%
  select(Cell,ChannelGroup)

dataScaled <- scale(data)
#mrr.dist <- dist(dataScaled, method = 'euclidean') #not needed as I use hclust,vector
#mrr.hclust.fit <- hclust(mrr.dist, method = 'complete') #causes memory overflow
tic("hclust")
mrrHclust.fit <- hclust.vector(dataScaled,method="ward",metric="euclidean")
toc()
mrrHclust.fit$notes <- "V4: Removed band & traffic. GSM900 Only"
saveRDS(mrrHclust.fit,file=paste0(exportPath,"HClustFitReducedFeats_V4_Nov26.rds"))

mrrHclust.fit <- readRDS(file=paste0(exportPath,"HClustFitReducedFeats_V4_Nov26.rds"))

# Dendogram plot
k=2
h=150
dendMrr <- as.dendrogram(mrrHclust.fit)
dendColored <- color_branches(dendMrr, h = h,groupLabels = TRUE)
plot(dendColored,leaflab="none")

# Cut the tree
mrr.clusters <- as.factor(cutree(mrrHclust.fit, h = h))
table(mrr.clusters)
mrrClusters.df <- bind_cols(labels,data.frame(cluster=mrr.clusters))
gmrrClustered <- inner_join(gmrr,mrrClusters.df,by=c("Cell", "ChannelGroup"))
table(gmrrClustered$cluster,gmrrClustered$Band) # ==> cluster 1 is mostly 1800, cluster is mostly 900



# Plot the clusters -------------------------------------------------------
# For plotting
mrrHclust.fit$notes
gmrrTidy <-  gmrrClustered %>%
  mutate_at(vars(`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`),funs(. / NoReportsPassedFilter)) %>%
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

map(kpiListBi1,plotBiDirLine)
map(kpiListBi2,plotBiDirBar)
map(kpiListUni1,plotUniDirLine)
map(kpiListUni1,plotUniDirLineC)
map(kpiListUni2,plotUniDirBar,12)


# Geographical Plots ------------------------------------------------------
atollCells <- readOGR(geoPath,"atollCells")
atollCells <- inner_join(atollCells,mrrClusters.df,by=c("CellName"="Cell"))
tmap_mode("view")
#myPal <- c("white","grey","grey","grey","green","grey","blue","red","grey")
myPal <- brewer.pal(8,"Set1")
tm_shape(atollCells) +
  tm_polygons("cluster",palette=myPal) 
