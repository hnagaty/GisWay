# WMRR reading and clustering
# The dfs used in below are
# wmrr ==> raw wmrr files
# wmrrRatio ==> from raw, but expressed in ratios instead of counts
# wmrrTidy ==> wmrrRatio reshaped for plotting. Also includes the cluster
# wmrrFeatures ==> wmrr quantiles, used as the clustering features

# Load libraries ----------------------------------------------------------
source("00_globalVars.R")
library(Hmisc)
library(fastcluster) # for hierarichal clustering, faster tha baser hclust
library(dendextend) # for plotting coloured dendrogram
library(tictoc) # for measuring time


# Some Functions ----------------------------------------------------------
# Given a kpi name, it get the exact name in the binMap
# e.g. "DL BLER" ==> "Bler"
getKpiName <- function(kpi, kpiList) {
  kpiList[which(map_lgl(kpiList, ~grepl(strpKpiName(.),strpKpiName(kpi))))] 
}

# Get the bin range value. The kpi name must not be an exact of that in the binMap
getBinRange <- function(bin, kpiname) {
  r = as.character(binMap[[getKpiName(kpiname,colnames(binMap))]][binMap$Bin==bin])
  if (length(r)==0) {return (NA)}
  return (r)
}

# same as above, but kpiname must be exact name as in binMap
getBinRange_ <- function(bin, kpishortname) {
  r = as.character(binMap[[kpishortname]][binMap$Bin==bin])
  if (length(r)==0) {return (NA)}
  return (r)
}

# Plotting functions
# All plotting functions work on wmrrTidy
# So, this df should be present in the parent environment before calling any of the fuctions.
plotCellWmrr <- function(cellname, service = "Speech AMR NB 12.2", crop = TRUE, save = FALSE) {
  wmrrCell <- wmrrTidy %>%
    filter(CellName == cellname, Service == service) %>%
    group_by(CellName, Service, Quantity, KpiShortName) %>%
    arrange(BinN) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(BinN)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(BinN) %>% 
    filter(!((cum1 == 0 | cum2 == 0) & crop)) %>%
    select(-cum1, -cum2)
  p <- ggplot(wmrrCell, aes(x=as.factor(Range), y=Value)) +
    geom_col(position = "dodge", fill="red4") +
    facet_wrap(~Quantity,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    scale_x_discrete("Bin Value") +
    ggtitle(paste0("WCDMA MRR for Cell ",cellname),subtitle = "WMRR") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  if (save) {
    ggsave(paste0("CellWMRR_",cellname,".png"), p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath))
  }
}

plotClusterWmrr <- function(cluster, service = "Speech AMR NB 12.2", pversion ="0", pdate = NULL, pnotes = NULL, crop = TRUE, save = TRUE) {
  wmrrCluster <- wmrrTidy %>%
    filter(Cluster == cluster, Service == service) %>%
    group_by(Service, Quantity, Cluster, BinN, Range, KpiShortName) %>%
    summarise(Value=mean(Value)) %>%
    arrange(BinN) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(BinN)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(BinN) %>% 
    filter(!((cum1 == 0 | cum2 == 0) & crop)) %>%
    select(-cum1, -cum2)
  ptitle <- paste0("WCDMA MRR, Cluster: ",cluster)
  psubtitle <- paste0("Version: ", pversion, " -> ", pdate, "\n", pnotes)
  p <- ggplot(wmrrCluster, aes(x=Range, y=Value)) +
    geom_col(position = "dodge", fill="red4") +
    facet_wrap(~Quantity,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(ptitle, psubtitle)
  print(p)
  if (save) {
    pfilename <- paste0("WMRRCluster_", pversion, "_Cluster_", cluster,".png")
    ggsave(pfilename, p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath,pfilename))
  }
}

plotKpiWmrr <- function(kpi, service = "Speech AMR NB 12.2", pversion ="0", pdate = NULL, pnotes = NULL, crop = TRUE, save = TRUE) {
  kpiS <- strpKpiName(kpi)
  kpiU <- gsub("[ /]","",kpi)
  wmrrKpi <- wmrrTidy %>%
    filter(strpKpiName(Quantity) == kpiS, Service == service) %>%
    group_by(Service, Quantity, Cluster, BinN, Range, KpiShortName) %>%
    summarise(Value=mean(Value)) %>%
    arrange(BinN) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(BinN)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(BinN) %>% 
    filter(!((cum1 == 0 | cum2 == 0) & crop)) %>%
    select(-cum1, -cum2)
  ptitle <- kpi
  psubtitle <- paste0("Version: ", pversion, " -> ", pdate, "\n", pnotes)
  p <- ggplot(wmrrKpi, aes(x=Range, y=Value)) +
    geom_col(position = "dodge", fill="red4") +
    facet_wrap(~Cluster,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(ptitle, psubtitle)
  print(p)
  if (save) {
    pfilename <- paste0("WMRR_KPI_", pversion, "_", kpiU,".png")
    ggsave(pfilename, p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath,pfilename))
  }
}


# Readin & data munging -----------------------------------------------
wmrr <- readMrrFiles(paths$wmrrDir,getMrrFiles(paths$wmrrDir))
# it's normal to have so many warnings in the WMRR, because not all lines have the same length

#Combine all instances into one row
wmrr <- wmrr %>%
  group_by(CellName,CellUserLabel,Service,Quantity) %>%
  summarise_all(sum) %>%
  ungroup()

# remove DL power from wmrr, as it's missing in many cells
# I may choose to add it later
wmrr <- wmrr %>%
  filter(!startsWith(Quantity,"DLT"))

# filter rows with all zeros
wmrr <- wmrr %>% 
  mutate(binsSum = rowSums(select(.,starts_with("Bin")), na.rm = TRUE)) %>%
  filter (binsSum !=0) %>%
  select(-CellName) %>%
  rename(CellName = CellUserLabel)

# Summarize the wmrr into quantiles
qtls <- seq(0.15,0.9,0.15) # the quantiles that are used to summarize the mrr data
qtlsnames <- paste0("Q",qtls*100)

# below df is needed for clustering (the features)
wmrrFeatures <- wmrr %>% # I re-used some code from gsm mrr; the code for computing quantiles
  filter (binsSum > 500) %>% # filter out cells with less than 1000 observations
  mutate_all(funs(replace(., is.na(.), 0))) %>%
  unite(cnt, starts_with("Bin"), sep=",") %>%
  mutate(val = paste0(seq (0,37, 1),collapse=",")) %>% # this step is not fool proof, if no. of bins is not 37
  rowwise() %>%
  mutate(q=paste0(wtd.quantile(as.numeric(strsplit(val,",",fixed=TRUE)[[1]]),
                               weights=as.numeric(strsplit(cnt,",",fixed=TRUE)[[1]]),
                               probs=qtls),
                  collapse=",")) %>%
  ungroup %>%
  select(-val,-cnt) %>%
  separate(q,sep=",",into=qtlsnames,convert=TRUE) %>%
  gather(key=qtl,value=value,-(CellName:Quantity)) %>%
  unite(ftr,Service, Quantity,qtl,sep="-") %>%
  spread(key=ftr,value = value)

# below is needed for plotting
wmrrRatio <- wmrr %>%
  select(-(Bin18:Bin37)) %>% # I dropped those bins as I already dropped the DL power
  mutate_at(vars(matches("^Bin[0-9]")), funs(./binsSum)) %>%
  rename(noSamples = binsSum)

# clustering
# identifier & notes for the model
data <- wmrrFeatures[,2:ncol(wmrrFeatures)]
labels <- wmrrFeatures[1]
cClass <- "WMRR"
cNotes <- "4th trial on WMRR. Euclidean distance, centroid linkage"
cVersion <- "v0.1d"
cDate <- as.Date("2018-12-18")
cList <- list(class=cClass, version=cVersion, date=cDate, notes=cNotes)
tic("Clustering Fit")
hClustFit <- hClustMrr(data, labels, cList, k=3, plotDendogram = TRUE)
toc()
hclusters <- hClustFit$clusters

# Explore other cuts
k <- 2
h <- NULL
newClusters <- tryHClust(model = hClustFit$model, labeled = hClustFit$clusters,
                         h = h, k =k)
hclusters <- newClusters

# Merge clusters with data & munging
wmrrTidy <- wmrrRatio %>% # for plotting purposes
  inner_join(hclusters) %>%
  gather(BinN, Value,starts_with("Bin")) %>%
  separate(BinN, into=c("dummy", "BinN"), sep="Bin", convert = TRUE) %>%
  select(-dummy) %>%
  rename(Cluster = cluster) %>%
  mutate(BinN = as.factor(BinN))
#rowwise() %>%
#mutate(Range = getBinRange(Bin, Quantity)) %>% #Very very slow & unstable
#ungroup()

# The below is not a straighforrwad approach, but it's because I intially wanted to introduced
# this directlty in the wmrrTidy
binMapDf <- wmrrTidy %>% distinct(Quantity,BinN) %>%
  arrange(Quantity,BinN) %>%
  rowwise() %>%
  mutate(KpiShortName = getKpiName(Quantity,colnames(binMap))) %>%
  mutate(Range = getBinRange_(BinN, KpiShortName)) %>%
  ungroup()
binMapDf$Range <- factor(binMapDf$Range, levels = unique(binMapDf$Range, ordered = TRUE))
wmrrTidy <- wmrrTidy %>%
  left_join(binMapDf)

# Save all data
saveRDS(wmrrFeatures, paste0(paths$exportPath,"wmrrFeatures_20181217.rds"))
saveRDS(wmrrRatio, paste0(paths$exportPath,"wmrrRatio_20181217.rds"))
saveRDS(wmrrTidy, paste0(paths$exportPath,"wmrrTidy_20181217.rds"))


# Load data ---------------------------------------------------------------
hClustFit <- readRDS("D:/Optimisation/~InProgress/201806_GisFramework/export/WMRR_HClustFit_v0.01b.rds")
cList <- hClustFit$model$info
wmrrFeatures <- readRDS(paste0(paths$exportPath,"wmrrFeatures_20181217.rds"))
wmrrRatio <- readRDS(paste0(paths$exportPath,"wmrrRatio_20181217.rds"))
wmrrTidy <- readRDS(paste0(paths$exportPath,"wmrrTidy_20181217.rds"))


# Example plots
getKpiName("dl bler", colnames(binMap))
plotCellWmrr("G5627B", crop = TRUE, save = FALSE)
plotKpiWmrr("DL Cpich Ecno")
plotClusterWmrr(8, crop = TRUE, save = TRUE,
                pversion = cList$version, pdate = cList$date,pnotes = cList$notes)

# Make plots of clusters
kpiList <- distinct(wmrrTidy,Quantity)$Quantity
walk(kpiList, plotKpiWmrr, pversion = cList$version, pdate = cList$date, pnotes = cList$notes)
clusterList <- sort(distinct(wmrrTidy,Cluster)$Cluster)
walk(clusterList, plotClusterWmrr, save=TRUE,
     pversion = cList$version, pdate = cList$date,pnotes = cList$notes)
