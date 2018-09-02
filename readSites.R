library(readr)
library(dplyr)
library(sp)
library(tidyr)
library(rgeos)
library(rgdal)
library(spdplyr)
library(ggplot2)
library(maptools)
library(RANN)
library(raster)
#library(distances)

# Some needed variables
globalPath <- "~/Dropbox/Voda/GISWay/sitesData/"
exportPath <- "~/Dropbox/Voda/GISWay/export/"
localPath <- "~/data/geodata/Siradel_20151217_Country50m/DLU/"
atollFile <- "atollDb.txt"
cellFile <- "oneCellDb.txt"
geoPath <- "/home/hnagaty/Dropbox/Voda/GISWay/geoData"
egUtmCRS <- CRS("+init=epsg:32636")
wgs84CRS <- CRS("+init=epsg:4326")

# Function to deduce site properties from the cell name
getPhySite <- function(cellname) {
  if (substr(cellname,1,1)=="L") {
    s=strsplit(cellname,"_",fixed=TRUE)
    if (length(s[[1]])==3) {
      pType=s[[1]][1]
      if (pType=="L18") {
        Type="LTE1800"
      } else if (pType=="L21") {
        Type="LTE2100"
      } else {
        Type="Invalid"
      }
      pSite=s[[1]][2]
      if ((nchar(pSite)==5 | nchar(pSite)==6) & substr(pSite,1,1)=="L") {
        PhySite <- substring(pSite,2)
      } else {
        PhySite="Invalid"
      }
    }
    else {
      PhySite <- "Invalid"
      Type <- "Invalid"
    }
  } else {
    chrLoc <- regexpr("[0-9]",cellname)[1]
    chrPrefix <- substring(cellname,1,chrLoc-1)
    numSuffix <- substring(cellname,chrLoc)
    if (nchar(chrPrefix)>2 | nchar(numSuffix)!=5) {
      PhySite="Invalid"
      Type="Invalid"
    } else {
      prefix <- substr(chrPrefix,1,1)
      suffix <- substr(numSuffix,1,4)
      if (prefix=="C") {
        PhySite <- paste0(substring(chrPrefix,2),suffix)
        Type <- "DCS1800"
      } else if (prefix=="G") {
        PhySite <- paste0(substring(chrPrefix,2),suffix)
        Type <- "UMTS2100"
      } else if (prefix=="M") {
        PhySite <- paste0(substring(chrPrefix,2),suffix)
        Type <- "UMTS900"
      } else {
        if (nchar(chrPrefix)<2) {
          PhySite <- paste0(chrPrefix,suffix)
          Type <- "GSM900" 
        } else {
          PhySite <- "Invalid"
          Type <- "Invalid"
        }
      }
    }     
  }
  if (PhySite!="Invalid") {
    a <- substr(cellname,nchar(cellname),nchar(cellname))
    Sector <- as.character(cellMap[a,2])
	if (Type=="UMTS2100") {Carrier <- paste0("G",as.character(cellMap[a,3]))}
	else if (Type="UMTS900") {Carrier <- paste0("M",as.character(cellMap[a,3]))}
	else {Carrier <- "NA"}
	} else {
    Sector="Invalid"
	  Carrier="Invalid"
  }
  return(paste(PhySite,Type,Sector,Carrier,sep=","))
}



cellMap <- read.delim("sectorMapping.txt")
cellMap$Carrier<-as.factor(cellMap$Carrier)
rownames(cellMap)<-cellMap$Suffix

atollDb <- read_delim(paste0(globalPath,atollFile),"\t",
                      escape_double = FALSE, trim_ws = TRUE, skip = 1,
                      col_names = c("Site","Lat","Long"),
                      col_types = "cdd")
# Remove duplicate sites
atollDb = atollDb[order(atollDb$Site),]
atollDb = atollDb[!duplicated(atollDb$Site),]

cellDb <- read_delim(paste0(globalPath,cellFile),"\t",
                     escape_double = FALSE, trim_ws = TRUE,col_types = "cc")
cellDb$Vendor <- as.factor(cellDb$Vendor)

cellDb <- cellDb %>%
  rowwise %>%
  mutate(siteData=getPhySite(CellName)) %>%
  separate(siteData,into=c("Site","Type","Sector"),sep=",") %>%    
  filter(Site!="Invalid") %>%
  ungroup()

siteDb <- cellDb %>%
  group_by(Site,Vendor) %>%
  summarise(noCells=n(),noSecs=n_distinct(Sector))

# join sites that exist in both atollDb and in CellDb
sites <- atollDb %>%
  inner_join(siteDb)

# Convert sites into sp object
coordinates(sites) <- c("Long","Lat")
proj4string(sites) <- wgs84CRS
utmSites <-spTransform(sites, egUtmCRS) # transformation is needed to convert units to meters

# Read voda regions
vodaRegions <- readOGR(geoPath,"Regions")
vodaRegions <-spTransform(vodaRegions, egUtmCRS)
vodaRegions <- vodaRegions %>%
  dplyr::select(REGION) %>%
  rename(Region=REGION)
vodaSubRegions <- readOGR(geoPath,"SubRegions")
vodaSubRegions <-spTransform(vodaSubRegions, egUtmCRS)
vodaSubRegions <- vodaSubRegions %>%
  dplyr::select(REGION,SUB_REGION) %>% 
  rename(Region=REGION,SubRegion=SUB_REGION)

# Read the administrative divisions
gSheakhat <- readOGR(geoPath,"SheakhatCleaned")
gSheakhat <-spTransform(gSheakhat, egUtmCRS)
gSheakhat <- gSheakhat %>%
  dplyr::select(SHYK_ANAME,SHYK_ENAME,QISM_ENAME,GOV_ENAME) %>%
  rename(SheakhaAr=SHYK_ANAME,SheakhaEn=SHYK_ENAME,QismEn=QISM_ENAME,GovernorateEn=GOV_ENAME)

# Read the geodata
siradel <- raster(paste0(localPath,"EGYPT_2_DLU_50m.bil"),crs="+init=epsg:32636")
siradelLegend <- read_delim(paste0(localPath,"EGYPT_2_DLU_50m.mnu"),
                            " ",col_names = c("Code","Clutter"))

# extract clutter type of each site
utmSites <- raster::extract(siradel,utmSites,sp=TRUE)
utmSites <- utmSites %>% rename(Clutter=EGYPT_2_DLU_50m)
utmSites@data$Clutter <- as.factor(utmSites@data$Clutter)
clutterNames <- siradelLegend$Clutter[siradelLegend$Clutter!="High_Dense_Vegetation"]
levels(utmSites@data$Clutter) <- clutterNames


# Spatial join
joinedSitesA <- over(utmSites,vodaSubRegions)
joinedSitesB <- over(utmSites,gSheakhat)
# Method #1
# the output is a data.frame
#sitesDf <- bind_cols(utmSites@data,joinedSitesA,joinedSitesB)
#write_csv(sitesDf,paste0(exportPath,"vodaSites.txt"))
# Method #2
# the output is an sp object
sitesSp <- spCbind(utmSites,joinedSitesA)
sitesSp <- spCbind(sitesSp,joinedSitesB)
sitesDf <- bind_cols(sitesSp@data,as.data.frame(sitesSp@coords))

# Next lines are for clustering the sites
sitesCord <- sitesDf %>%
  select(Long,Lat)
sitesNames <- sitesDf$Site
k=6
noSites <- NROW(sitesNames)
# not used library
#sitesDistances <- distances(sitesDf,id_variable = "Site",dist_variables = c("Long","Lat"))
#nearestSitesIdx <- nearest_neighbor_search(sitesDistances,k=k)
#nearestSitesNames <- matrix(sitesNames[nearestSitesIdx],nrow=k,ncol=noSites)

# K nearest sites
distMatrix <- nn2(sitesCord,k=k)
sitesDf$minDist <- apply(distMatrix[["nn.dists"]][,2:k],1,min)
sitesDf$maxDist <- apply(distMatrix[["nn.dists"]][,2:k],1,max)
sitesDf$meanDist <- apply(distMatrix[["nn.dists"]][,2:k],1,mean)
sitesDf$sdDist <- apply(distMatrix[["nn.dists"]][,2:k],1,sd)
sitesDf <- sitesDf %>% mutate(normSd=sdDist/meanDist)

ggplot(data=sitesDf,aes(x=Long,y=Lat,col=Vendor)) + geom_point()
qplot(Long,Lat,data=sitesDf,colour=Region)

# Clustering
# I still need to try to do clutering based on underlying clutter class
#sitesCl <- sitesDf %>% dplyr::select(Clutter,meanDist,sdDist)

#write_csv(sites,"vodaSitesDf.txt")

# Clustering
sitesCl <- sitesDf %>%
  dplyr::select(meanDist,normSd)

# Hierarchical Clustering
hclust.dist <- dist(sitesCl, method = "euclidean") # distance matrix
hclust.fit <- hclust(hclust.dist, method="ward.D2")
# For loop for different k values
sites.clusters <- data.frame(Site=character(),
                             Long=double(),
                             Lat=double(),
                             Cluster=integer(),
                             kValue=integer(),
                             stringsAsFactors=FALSE)
sitesLoc <- dplyr::select(sitesDf,Site,Long,Lat)
for (k in 2:7) {
  clusters <- cutree(hclust.fit, k=k)
  sitesLoc$Cluster <- clusters
  sitesLoc$kValue <- k
  sites.clusters <- bind_rows(sites.clusters,sitesLoc)
}
sites.clusters$kValue <- as.factor(sites.clusters$kValue)
sites.clusters$Cluster <- as.factor(sites.clusters$Cluster)
p <- ggplot(data=sites.clusters,aes(x=Long,y=Lat)) +
  facet_wrap(~kValue) +
  geom_point()
p + aes(col=Cluster)

for (k in 2:7) {
  sites.export <- dplyr::filter(sites.clusters,kValue==k)
  write_csv(sites.export,paste0("clusteredSitesK",k,".txt"))
}
write_csv(sites.clusters,"clusteredSitesAllK.txt")

# Using k=5
sites.5clusters <- dplyr::filter(sites.clusters,kValue==5)
levels(sites.5clusters$Cluster) <- c("Urban","Rural","Road","Isolated","FarIsolated","NA","NA")
sites.5clusters <- sites.5clusters %>%
  select(Site,Cluster)
# Merge with the sitesDf
sitesDf <- inner_join(sitesDf,sites.5clusters)
write_csv(sitesDf,paste0(exportPath,"vodaSites.txt"))
table(sitesDf$Region,sitesDf$Cluster)
