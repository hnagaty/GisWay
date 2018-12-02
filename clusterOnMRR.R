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
library(rlist)

# Parallel Processing -----------------------------------------------------
registerDoParallel(cores=2) #multi-core like functionality


# Functions & Indentifiers ------------------------------------------------
pasteDir <- function(c) { # needed for GMRR parsing
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}

myDesktop <- "Linux"
#myDesktop <- "Windows"

if (myDesktop == "Linux") {
  dataDir <- "~/data/gmrr/"
  exportPath <- "~/Dropbox/Voda/GISWay/export/"
  pngPath <- "~/Dropbox/Voda/GISWay/gisWay/outs/plots/"  
} else {
  dataDir<- "d:/data/mrr/2018Sep25/"
  exportPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
  pngPath <- "D:/Optimisation/~InProgress/201806_GisFramework/gisWay/outputs/plots/"
  geoPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
}


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
gmrrRatio <- gmrr %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter)) %>%
  select(-contains("TrafficLevelM"))

write_csv(gmrr,paste0(exportPath,"GsmMRR.csv"))
write_csv(gmrrRatio,paste0(exportPath,"GsmMrrRatio.csv"))

# For Clustering II -----------------------------------------------------------
# read data
gmrr <- read_csv(paste0(exportPath,"GsmMRR.csv"))
gmrrRatio <- read_csv(paste0(exportPath,"GsmMrrRatio.csv"))
nCols <- ncol(gmrr)

# define stops
qtls <- seq(0.15,0.9,0.15)
qtlsnames <- paste0("Q",qtls*100)

# reduce features
tic("GMRR Features")
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

gmrrFeatures <- readRDS(file=paste0(exportPath,"gmrrFeatures.rds"))

# data & labels
#excludedFeatures <- vars(Band,traffic) # I couldn't use this, but I want to.
data <- gmrrFeatures %>%
  filter(Band=="GSM900") %>%
  select(-Cell,-ChannelGroup,-Band,-traffic)
labels <- gmrrFeatures %>%
  filter(Band=="GSM900") %>%
  select(Cell,ChannelGroup)

cNotes <- "Some notes in here"
cVersion <- "v31"
cDate <- "2018-12-03"
cList <- list(version=cVersion, date=cDate, notes=cNotes)

# clustering
hClustMrr <- function (data, infoList, save=TRUE, plotDendogram=TRUE) {
  k=8
  h=180
  scaled <- scale(data)
  fit <- hclust.vector(scaled,method="ward",metric="euclidean")
  fit <- list.append(fit,info=infoList)
  if (save) saveRDS(fit,file=paste0(exportPath,"HClustFit_",infoList$version,".rds"))
  # Dendogram plot
  if (plotDendogram) {
    dendMrr <- as.dendrogram(fit)
    dendColored <- color_branches(dendMrr, h = h) #,groupLabels = TRUE) ==> seems to not match clusters in cuttree
    plot(dendColored,leaflab="none")
  }
  # Cut the tree
  clusters <- as.factor(cutree(mrrHclust.fit, k = k))
  return(list(model=fit,clusters=clusters))
}

cd1 <- hClustMrr(data,cList)

mrrClusters.df <- bind_cols(labels,data.frame(cluster=mrr.clusters))
gmrrClustered <- inner_join(gmrr,mrrClusters.df,by=c("Cell", "ChannelGroup"))
table(gmrrClustered$cluster,gmrrClustered$Band) #


mrrHclust.fit <- readRDS(file=paste0(exportPath,"HClustFit_Nov21.rds"))


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
myPal <- c("grey","grey","blue","grey","grey","grey","grey","grey")
#myPal <- brewer.pal(8,"Set1")
tm_shape(atollCells) +
  tm_polygons("cluster",palette=myPal) 

gmrrTidy %>% filter(cluster==5) %>%
  ggplot(aes(x=Bin,y=Value,fill=dir)) +
  geom_col(position = "dodge",col="grey") +
  facet_wrap(~kpi,ncol=3,scales="free_x", space="free_x") +
  scale_y_continuous(labels = scales::percent)

gmrrClustered %>% filter(Cell=="58874") %>% select(cluster)

# Single cell plot

# Plot Cell MRR Function
  # plot mrr for a single cell
  # the gmrrRatio df should be present
plotCellMrr <- function(neededCell, save=FALSE) {
  gmrrCell <-  gmrrRatio %>%
    filter(Cell==neededCell) %>%
    gather("Measure","Value",`RXQUALUL(0,0)`:`PATHLOSSDIFF(25,25)`) %>%
    select(BSC:Band,Measure,Value) %>%
    separate(Measure,into=c("Measure","Bin","c","d"),sep="[\\)\\(,]",convert=TRUE) %>%
    select(-c,-d) %>%
    mutate(dir=substring(Measure,regexpr("UL|DL",Measure))) %>%
    separate(Measure,into=c("kpi"),sep="UL|DL") %>%
    rowwise() %>%
    mutate(dir = pasteDir(dir)) %>%
    ungroup() %>%
    select(BSC:kpi, dir, Bin, Value, -SubCellType, -Band, -BSC) %>%
    mutate(Bin=ifelse(kpi=="RXLEV",Bin-110,Bin)) %>%
    group_by(Cell, ChannelGroup, kpi, dir) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(Bin)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(Bin) %>% 
    filter(cum1 != 0 & cum2 != 0)
  p <- ggplot(gmrrCell, aes(x=Bin, y=Value,fill=dir)) +
    geom_col(position = "dodge") +
    facet_wrap(~kpi,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    ggtitle(paste0("GSM MRR for cell ",neededCell),subtitle = "notes go in here")
  if (save) {
    ggsave(paste0("V4_Cell_",neededCell, ".png"), p, path=pngPath,width=30, height=15, units="cm")    
  }
  p
}

cells <- list("55061","55063","D31911")
map(cells, plotCellMrr, save=TRUE)

