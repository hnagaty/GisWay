# WMRR reading and clustering

# Load libraries ----------------------------------------------------------
source("00_globalVars.R")
library(Hmisc)
library(fastcluster) # for hierarichal clustering, faster tha baser hclust
library(dendextend) # for plotting coloured dendrogram
library(tictoc) # for measuring time


# Read in the WCDMA MRR files -----------------------------------------------
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

#wmrrTest <- wmrr %>%
#  separate(CellName, into = c("RNC", "C"), sep = "/") %>%
#  select(-C) %>%
#  rename(CellName = CellUserLabel) %>%
#  filter(CellName == "G43481") %>%
#  gather(BinC, Value, starts_with("Bin")) %>%
#  filter (!is.na(Value)) %>%
#  mutate(Bin = as.integer(sub("Bin","", BinC))) %>%
#  select(-BinC) %>%
#  arrange(Quantity,Bin)

# Summarize the wmrr into quantiles
qtls <- seq(0.15,0.9,0.15) # the quantiles that are used to summarize the mrr data
qtlsnames <- paste0("Q",qtls*100)

wmrrFeatures <- wmrr %>% # I re-used some code from gsm mrr; the code for computing quantiles
  #filter (CellName=="GD3153A") %>%
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


# clustering
# identifier & notes for the model
data <- wmrrFeatures[,2:ncol(wmrrFeatures)]
labels <- wmrrFeatures[1]
cClass <- "WMRR"
cNotes <- "1st b trial on WMRR"
cVersion <- "v0.01b"
cDate <- as.Date("2018-12-16")
cList <- list(class=cClass, version=cVersion, date=cDate, notes=cNotes)
tic("Clustering Fit")
hClustFit <- hClustMrr(data, labels, cList, k=8)
toc()
hClustFit <- readRDS("D:/Optimisation/~InProgress/201806_GisFramework/export/WMRR_HClustFit_v0.01b.rds")
cList <- hClustFit$model$info

# Explore other cuts
k <- NULL
h <- 400
newClusters <- tryHClust(model = hClustFit$model, labeled = hClustFit$clusters,
                         h = h, k =k)
wmrrRatio <- wmrr %>%
  select(-(Bin18:Bin37)) %>% # I dropped those bins as I already dropped the DL power
  mutate_at(vars(matches("^Bin[0-9]")), funs(./binsSum)) %>%
  rename(noSamples = binsSum)

hclusters <- hClustFit$clusters #newClusters
wmrrTidy <- wmrrRatio %>% # for plotting purposes
  inner_join(hclusters) %>%
  gather(Bin, Value,starts_with("Bin")) %>%
  separate(Bin, into=c("dummy", "Bin"), sep="Bin", convert = TRUE) %>%
  select(-dummy) %>%
  rename(Cluster = cluster)

plotCellWmrr <- function(cellname, service = "Speech AMR NB 12.2", crop = TRUE, save = FALSE) {
  wmrrCell <- wmrrTidy %>%
    filter(CellName == cellname, Service == service) %>%
    group_by(CellName, Service, Quantity) %>%
    arrange(Bin) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(Bin)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(Bin) %>% 
    filter(!((cum1 == 0 | cum2 == 0) & crop)) %>%
    select(-cum1, -cum2)
  
  p <- ggplot(wmrrCell, aes(x=Bin, y=Value)) +
    geom_col(position = "dodge", fill="red4") +
    facet_wrap(~Quantity,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    ggtitle(paste0("WCDMA MRR for Cell ",cellname),subtitle = "WMRR")
  print(p)
  if (save) {
    ggsave(paste0("CellWMRR_",cellname,".png"), p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath))
  }
}

plotCellWmrr("G43481", crop = FALSE, save = TRUE)


plotClusterWmrr <- function(cluster, service = "Speech AMR NB 12.2", pversion ="0", pdate = NULL, pnotes = NULL, crop = TRUE, save = TRUE) {
  wmrrCluster <- wmrrTidy %>%
    filter(Cluster == cluster, Service == service) %>%
    group_by(Service, Quantity, Cluster, Bin) %>%
    summarise(Value=mean(Value)) %>%
    arrange(Bin) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(Bin)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(Bin) %>% 
    filter(!((cum1 == 0 | cum2 == 0) & crop)) %>%
    select(-cum1, -cum2)
  ptitle <- paste0("WCDMA MRR, Cluster: ",cluster)
  psubtitle <- paste0("Version: ", pversion, " -> ", pdate, "\n", pnotes)
  p <- ggplot(wmrrCluster, aes(x=Bin, y=Value)) +
    geom_col(position = "dodge", fill="red4") +
    facet_wrap(~Quantity,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    ggtitle(ptitle, psubtitle)
  print(p)
  if (save) {
    pfilename <- paste0("WMRRCluster_", pversion, "_Cluster_", cluster,".png")
    ggsave(pfilename, p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath,pfilename))
  }
}

plotClusterWmrr(8, crop = TRUE, save = TRUE,
               pversion = cList$version, pdate = cList$date,pnotes = cList$notes)

clusterList <- sort(distinct(wmrrTidy,Cluster)$Cluster)

walk(clusterList, plotClusterWmrr, save=TRUE,
    pversion = cList$version, pdate = cList$date,pnotes = cList$notes)




plotKpiWmrr <- function(kpi, service = "Speech AMR NB 12.2", pversion ="0", pdate = NULL, pnotes = NULL, crop = TRUE, save = TRUE) {
  kpiS <- strpKpiName(kpi)
  kpiU <- gsub("[ /]","",kpi)
  wmrrKpi <- wmrrTidy %>%
    filter(strpKpiName(Quantity) == kpiS, Service == service) %>%
    group_by(Service, Quantity, Cluster, Bin) %>%
    summarise(Value=mean(Value)) %>%
    arrange(Bin) %>%
    mutate(cum1=cummax(Value)) %>%
    arrange(desc(Bin)) %>%
    mutate(cum2=cummax(Value)) %>%
    arrange(Bin) %>% 
    filter(!((cum1 == 0 | cum2 == 0) & crop)) %>%
    select(-cum1, -cum2)
  ptitle <- kpi
  psubtitle <- paste0("Version: ", pversion, " -> ", pdate, "\n", pnotes)
  p <- ggplot(wmrrKpi, aes(x=Bin, y=Value)) +
    geom_col(position = "dodge", fill="red4") +
    facet_wrap(~Cluster,ncol=3,scales="free") +
    scale_y_continuous(name=NULL, labels = scales::percent) +
    ggtitle(ptitle, psubtitle)
  print(p)
  if (save) {
    pfilename <- paste0("WMRR_KPI_", pversion, "_", kpiU,".png")
    ggsave(pfilename, p, path=paths$pngPath,width=30, height=15, units="cm")
    message(paste0("Image saved in ",paths$pngPath,pfilename))
  }
}


plotKpiWmrr("DL  Cpich Ecno")

kpiList <- distinct(wmrrTidy,Quantity)$Quantity

walk(kpiList, plotKpiWmrr)

#===========================================================================================
