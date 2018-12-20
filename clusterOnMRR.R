# I should run the file "readGenericMRR.R" before this
# the above named script reads the names of complete MRR export files in a given folder

source("00_globalVars.R")

# Load libraries ----------------------------------------------------------
library(Hmisc) # for computation of the quantiles
library(fastcluster) # for hierarichal clustering, faster tha baser hclust
library(dendextend) # for plotting coloured dendrogram
library(tictoc) # for measuring time
library(dummies)
library(sp)
library(rgdal)
library(spdplyr)
library(tmap)
library(RColorBrewer)
library(viridisLite)
library(hablar) # for I don't know what
library(readxl)

# Environment setup -------------------------------------------------------
#registerDoParallel(cores=2) #multi-core like functionality

# Load some data

# Functions & Indentifiers ------------------------------------------------
# set init value in parsing of CNA GSM parameters
setInitValue <- function(x,param) {
  x[is.na(x)] <- pull(defaultParams[defaultParams[1]==param,2])
  retype(x)
}


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
    ggtitle(paste0("GSM MRR for cell ",neededCell),subtitle = "After")
  if (save) {
    ggsave(paste0("CellMRR_",neededCell, "_02After.png"), p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath))
  }
  p
}


#Plotting functions
kpiListBi1 <- c("PATHLOSS","RXLEV")
kpiListBi2 <- c("RXQUAL")
kpiListUni1 <- c("PATHLOSSDIFF")
kpiListUni2 <- c("BSPOWER","MSPOWER","TAVAL")

#Generic plotting functions
plotBiDirLine <- function(df, kpiV="RXLEV", subtitle=NULL) {
    filter(df, kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=dir)) +
    geom_line(position = "dodge",size=1) +
    facet_wrap(~cluster,ncol=3) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle=subtitle,x="Value",y="Percentage",color="Direction")
}

plotBiDirBar <- function(df, kpiV="RXQUAL",  subtitle=NULL) {
    filter(df, kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,fill=dir)) +
    geom_col(position = "dodge",col="grey") +
    facet_wrap(~cluster,ncol=3) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle=subtitle,x="Value",y="Percentage",fill="Direction")
}

plotUniDirLine <- function(df, kpiV="PATHLOSSDIFF",  subtitle=NULL) {
    filter(df, kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=cluster)) +
    geom_line(position = "identity",size=1) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle=subtitle,x="Value",y="Percentage",color="Direction")
}

plotUniDirLineC <- function(df, kpiV="PATHLOSSDIFF",  subtitle=NULL) {
    filter(df, kpi==kpiV) %>%
    ggplot(aes(x=Bin,y=Value,col=cluster)) +
    facet_wrap(~cluster,ncol=3) +
    geom_line(position = "identity",size=1) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle=subtitle,x="Value",y="Percentage",color="Direction")
}

plotUniDirBar <- function(df, kpiV="TAVAL",taLimit=20,  subtitle=NULL) {
    filter(df, kpi==kpiV) %>%
    filter(!(kpi=="TAVAL" & Bin >= taLimit )) %>%
    #mutate(Bin=as.factor(Bin)) %>%
    ggplot(aes(x=Bin,y=Value)) +
    geom_col(fill="orange") +
    facet_wrap(~cluster,ncol=3) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=kpiV,subtitle=subtitle,x="Value",y="Percentage",fill="Direction")
}

plotCluster <- function(df, cl, subtitle=NULL) {
  filter(df, cluster==cl) %>%
    ungroup() %>%
    select(-cluster) %>%
    group_by(kpi, dir) %>%
    arrange(desc(Bin)) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(Bin) %>% 
    mutate(cum2=cummax(Value)) %>%
    filter(cum1 > 1e-03 & cum2 > 1e-03) %>%
  ggplot(aes(x=Bin, y=Value,fill=dir)) +
    geom_col(position = "dodge") +
    facet_wrap(~kpi,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) + 
    labs(title=paste0("Cluster ",cl),subtitle=subtitle,x=NULL,y="Percentage",fill="Direction")
}

# Plot the clusters
plotMrrClusters <- function(df, save=FALSE, filename=NULL, subtitle=NULL,byCluster=FALSE) {
  cnts <- count(gmrrClustered,cluster)
  dfP <-  df %>%
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
    mutate(Bin=ifelse(kpi=="RXLEV",Bin-110,Bin)) %>%
    group_by(cluster, kpi, dir) %>%
    arrange(desc(Bin)) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(Bin) %>% 
    mutate(cum2=cummax(Value)) %>%
    filter(cum1 > 1e-05 & cum2 > 1e-05)
  if (byCluster) { # each cluster in a separate page
    noClusters <- n_distinct(df$cluster)
    p <- map2((rep(list(dfP),noClusters)), as.list(1:noClusters), plotCluster, subtitle=subtitle)
    if (save) {
      imap(p,~ggsave(filename=paste0(filename,"Cluster",.y,".png"),
                     plot=., path=pngPath, width=30, height=15, units="cm"))
    }
    p
  } else { # each kpi in a separate page
    p1 <- map2(rep(list(dfP),length(kpiListBi1)),  kpiListBi1, plotBiDirLine, subtitle=subtitle)
    p2 <- map2(rep(list(dfP),length(kpiListBi2)),  kpiListBi2, plotBiDirBar, subtitle=subtitle)
    p3 <- map2(rep(list(dfP),length(kpiListUni1)), kpiListUni1, plotUniDirLine, subtitle=subtitle)
    p4 <- map2(rep(list(dfP),length(kpiListUni1)), kpiListUni1, plotUniDirLineC, subtitle=subtitle)
    p5 <- map2(rep(list(dfP),length(kpiListUni2)), kpiListUni2, plotUniDirBar, subtitle=subtitle)
    ps <- list(p1, p2, p3, p4, p5)
    if (save) {
        map(flatten(ps),~ggsave(filename=paste0(filename,"KPI_",.$labels$title,".png"),
                                plot=., path=pngPath, width=30, height=15, units="cm"))
    }
    ps
  }
}

# Read MRR files ----------------------------------------------------------
# mrrconf.df should be pre-defined. It's defined in "readGenericMRR.R"
gmrrDir <- "D:/Optimisation/~InProgress/201806_GisFramework/ossData/mrr_cluster7_after/"
gmrr <- readMrrFiles(gmrrDir,getMrrFiles(gmrrDir))

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

write_csv(gmrr,paste0(paths$exportPath,"GsmMRR.csv"))
gmrr <- read_csv(paste0(paths$exportPath,"GsmMRR.csv"))

nCols <- ncol(gmrr)
gmrrRatio <- gmrr %>%
  mutate_at(11:nCols,funs(. / NoReportsPassedFilter)) %>%
  select(-contains("TrafficLevelM"))

# reduce features
qtls <- seq(0.15,0.9,0.15) # the stops
qtlsnames <- paste0("Q",qtls*100)

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

# clustering
# identifier & notes for the model
cNotes <- "900 band only. 8 clusters. Aggregated quantiles (6 quantiles).This is the version where subsequent analysis is made"
cVersion <- "v4"
cDate <- as.Date("2018-12-09")
cList <- list(version=cVersion, date=cDate, notes=cNotes)
tic("Clustering Fit")
hClustFit <- hClustMrr(data, labels, cList, k=8)
toc()

# Explore other cuts
k <- 8
h <- NULL
newClusters <- tryHClust(model = hClustFit$model, labeled = hClustFit$clusters,
                         h = h, k =k)
# Now, merge the final clusters into the original dataframe
# Note that there is potential mismatch in the saved files here,
# as the txt/csv data is saved within hClustFit, and the images are saved based on clusters of tryHclust
rm(gmrrClustered)
gmrrClustered <- inner_join(gmrr,newClusters,by=c("Cell", "ChannelGroup"))
table(gmrrClustered$cluster,gmrrClustered$Band)


plotMrrClusters(gmrrClustered, save=TRUE, byCluster = TRUE,
                filename=paste(cList$version,cList$date,"", sep="_"),
                subtitle=paste(cList$version, cList$notes, sep="\n"))


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


cells <- list("55061","55063","D31911")
map(cells, plotCellMrr, save=TRUE)



# Read other supplementary data -------------------------------------------
vfSites.df <- read_csv(paste0(exportPath,"vodaSites_201811b.csv"))

gsmKpi <- read_tsv(paste0(ossDataDir,"gsmKpi_122018.txt"),
                   skip=1, na=c("","NA","#DIV/0"))

gsmKpi2 <- gsmKpi %>%
  mutate(IcmAll = IcmB1+IcmB2+IcmB3+IcmB4+IcmB5) %>%
  mutate_at(vars(starts_with("IcmB")), funs(./IcmAll))

pwrCntrlParams <- read_xlsx(paste0(ossDataDir,"cnadb/pwrControl_20181204.xlsx"),
                             col_types = c("text", "text", "numeric","numeric", "numeric",
                                           "numeric", "numeric", "numeric", "numeric",
                                           "numeric", "numeric", "numeric", "numeric",
                                           "numeric", "numeric", "numeric", "numeric",
                                           "numeric", "numeric", "numeric", "numeric",
                                           "numeric", "numeric", "text", "numeric"))

defaultParams <- read_xlsx(paste0(ossDataDir, "cnadb/defaultParams.xlsx"))

for (i in seq(3,ncol(pwrCntrlParams))) {
  pwrCntrlParams[i] <- setInitValue(pwrCntrlParams[[i]],names(pwrCntrlParams)[i])
}

remedyIncs <- read_csv(paste0(miscDataDir,"remedyIncs_20181202.csv"),
                       col_names = c("IncID", "AssetID", "Status", "FaultTime", 
                                     "ResTime", "ResFault", "Desc"))
sitesFaults <- count(remedyIncs,AssetID)
hccClusters <- hClustFit$clusters %>%
  rowwise %>%
  mutate(siteData=getPhySite(Cell)) %>%
  separate(siteData,into=c("Site","Type","Sector","Carrier"),sep=",")
clustersIncs <- inner_join(hccClusters,sitesFaults, by = c("Site" = "AssetID"))
ggplot(clustersIncs,aes(x=n)) +
  geom_histogram(binwidth = 5) +
  facet_wrap(~cluster, scale="free_y") +
  scale_x_continuous(limits = c(0,100))

clusterParam <- inner_join(pwrCntrlParams,hClustFit$clusters)
ggplot(clusterParam,aes(x=ssdesdl)) +
  geom_histogram(binwidth=1) +
  facet_wrap(~cluster, scale="free_y") 
  scale_x_continuous(breaks=seq(80,95,1))

clusterKpi <- inner_join(gsmKpi2, hClustFit$clusters, by = c("CellName" = "Cell"))
ggplot(clusterKpi,aes(x=cluster,y=IcmB1)) +
  geom_boxplot()


clusterSiteInfo <- vfSites.df %>%
  select(Site,Region, SubRegion, SiteClass, meanDist, maxDist) %>%
  right_join(hccClusters)

ggplot(clusterSiteInfo, aes(x=cluster, y=meanDist)) + geom_boxplot()
ggplot(clusterSiteInfo,aes(x=meanDist)) +
  geom_histogram() + 
  facet_wrap(~cluster, scales = "free_y")

table(clusterSiteInfo$SiteClass, clusterSiteInfo$cluster)

#====================================================

dfP <-  gmrrClustered %>%
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
  mutate(Bin=ifelse(kpi=="RXLEV",Bin-110,Bin)) %>%
  group_by(cluster, kpi, dir) %>%
  arrange(desc(Bin)) %>%
  mutate(cum1=cummax(Value)) %>%
  arrange(Bin) %>% 
  mutate(cum2=cummax(Value)) %>%
  filter(cum1 > 1e-05 & cum2 > 1e-05)

filename <- "Hany86"
p <- map2((rep(list(dfP),noClusters)), as.list(1:noClusters), plotCluster)



plotCellMrr("D80301", save = TRUE)

mrrClusters <- read_csv(paste0(paths$exportPath, "HClustOutput_v4_clusters.csv"))

clust7 <- filter(mrrClusters, cluster == 7) %>% distinct(Cell)
clust7l <- clust7[[1]]

clust7 <- read_csv("cluster7indelta.csv")
clust7l <- clust7[[1]]

map(clust7l, safely(plotCellMrr), save = TRUE)

