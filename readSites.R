# Remebering the variable names
# atollDb --> The Atoll DB (per PhySite, live only)
# cellDb --> The OneCell Database
# siteDb --> The site database, as extracted from the OneCell DB
# vfSites --> The result after joining OneCell DB & Atoll DB. This is a projected (UTM Zone 36) SpatialPointDataFrame
#             This is the main data structure
# vfSites.df --> df of vfSites, with the coords in lat/long

# Load Libraries ----------------------------------------------------------
library(raster)
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
library(tmap)
library(ggmap)
library(dendextend)

# Define variables & functions --------------------------------------------
# Variables
#parentPath <- "~/Dropbox/Voda/GISWay/"
parentPath <- "D:/Optimisation/~InProgress/201806_GisFramework/"
mainPath <- paste0(parentPath,"gisWay/")
globalPath <- paste0(parentPath,"sitesData/")
exportPath <- paste0(parentPath,"export/")
geoPath <- paste0(parentPath,"geoData")
#localPath <- "~/data/geodata/Siradel_2016March/DLU/"
localPath <- "D:/Optimisation/~InProgress/201806_GisFramework/geoData/Siradel_2016_Country/DLU/"
atollFile <- "atollSites.csv"
cellFile <- "oneCell.csv"

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
      } else if (pType=="L09") {
        Type="LTE900"}
      else {
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
	if (Type=="UMTS2100" & !is.na(Sector)) {Carrier <- paste0("G",as.character(cellMap[a,3]))}
	else if (Type=="UMTS900" & !is.na(Sector)) {Carrier <- paste0("M",as.character(cellMap[a,3]))}
	else {Carrier <- "NA"}
	} else {
    Sector="Invalid"
	  Carrier="Invalid"
  }
  return(paste(PhySite,Type,Sector,Carrier,sep=","))
}



# Read in the data --------------------------------------------------------
# Cellmap file, maps cell suffix to sector & carrier. This is not the latest version
cellMap <- read.delim(paste0(mainPath,"sectorMapping.txt"))
cellMap$Carrier<-as.factor(cellMap$Carrier)
rownames(cellMap)<-cellMap$Suffix

atollDb <- read_csv(paste0(globalPath,atollFile),
                    trim_ws = TRUE, skip = 1,
                    col_types=cols(
                      SiteID = col_skip(),
                      PhySite = col_character(),
                      NAME = col_skip(),
                      LATITUDE = col_double(),
                      LONGITUDE = col_double(),
                      NAME_1 = col_skip(),
                      SITE_TYPE = col_skip()
                    ))
colnames(atollDb) <- c("Site","Lat","Long")
# I didn't investigate the warnings. Don't know why they are here
atollDb <- atollDb %>%
  filter(complete.cases(.)) %>%
  distinct(Site,.keep_all=TRUE)

# Remove duplicate sites, the R way
#atollDb = atollDb[order(atollDb$Site),]
#atollDb = atollDb[!duplicated(atollDb$Site),]

cellDb <- read_csv(paste0(globalPath,cellFile),
                   trim_ws = TRUE,
                   col_types = cols_only(
                     CELL = col_character(),
                     Vendor = col_factor(levels=c("ER","HU"))
                   ))
colnames(cellDb) <- c("CellName","Vendor")

# Deduce site info from the cell name sufix
cellDb <- cellDb %>%
  rowwise %>%
  mutate(siteData=getPhySite(CellName)) %>%
  separate(siteData,into=c("Site","Type","Sector","Carrier"),sep=",") %>%    
  filter(Site!="Invalid") %>%
  ungroup()

#test <- cellDb %>%  filter(Site=="2999")

# the below returns duplicate sites if a site doesn't have same no. of sectors for each technology
# it needs revisiting -- DONE
siteDb <- cellDb %>%
  group_by(Site,Vendor,Type) %>%
  summarise(noCells=n(),noSecs=n_distinct(Sector)) %>%
  spread(Type,noCells,fill=0,sep="_noCells") %>%
  group_by(Site,Vendor) %>%
  summarise_all(max)
colnames(siteDb) <- gsub("Type_","",colnames(siteDb))
# there is one duplicate site (site 2999), that is Hu & Er
# I didn't take an action into it
# Should classify it as "Dual Vendor"

# join sites that exist in both atollDb and in CellDb
vfSites <- atollDb %>%
  inner_join(siteDb) %>%
  arrange(Vendor,Site)

# Convert sites into sp object
coordinates(vfSites) <- c("Long","Lat")
proj4string(vfSites) <- wgs84CRS
vfSites <-spTransform(vfSites, egUtmCRS) # projecting to utm zone 36 (convert units to meters)


# Read in the Geodata -----------------------------------------------------
# Read voda regions
vfRegions <- readOGR(geoPath,"Regions")
vfRegions <-spTransform(vfRegions, egUtmCRS)
vfRegions <- vfRegions %>%
  dplyr::select(REGION) %>%
  rename(Region=REGION)
vfSubRegions <- readOGR(geoPath,"SubRegionss")
vfSubRegions <-spTransform(vfSubRegions, egUtmCRS)
vfSubRegions <- vfSubRegions %>%
  dplyr::select(REGION,SUB_REGION) %>% 
  rename(Region=REGION,SubRegion=SUB_REGION)

# Read the administrative divisions
#gSheakhat <- readOGR(geoPath,"SheakhatCleaned")
gSheakhat <- readOGR(geoPath,"SHYK_2014Jan")
gSheakhat <-spTransform(gSheakhat, egUtmCRS)
gSheakhat <- gSheakhat %>%
  dplyr::select(SHYK_ANAME,SHYK_ENAME,QISM_ENAME,GOV_ENAME) %>%
  rename(SheakhaAr=SHYK_ANAME,SheakhaEn=SHYK_ENAME,QismEn=QISM_ENAME,GovernorateEn=GOV_ENAME)

# Read the geodata
siradel <- raster(paste0(localPath,"EGYPT_3_DLU_50m.bil"),crs="+init=epsg:32636")
siradelLegend <- read_delim(paste0(localPath,"EGYPT_3_DLU_50m.mnu"),
                            " ",col_names = c("Code","Clutter"))
# extract clutter type of each site
vfSites <- raster::extract(siradel,vfSites,sp=TRUE)
vfSites <- vfSites %>% rename(Clutter=EGYPT_3_DLU_50m)
vfSites@data$Clutter <- as.factor(vfSites@data$Clutter)
clutterNames <- siradelLegend$Clutter[siradelLegend$Clutter!="High_Dense_Vegetation"]
levels(vfSites@data$Clutter) <- clutterNames


# Spatial join
joinedSitesA <- over(vfSites,vfSubRegions)
joinedSitesB <- over(vfSites,gSheakhat)
# Method #1
# the output is a data.frame
#sitesDf <- bind_cols(vfSites@data,joinedSitesA,joinedSitesB)
#write_csv(sitesDf,paste0(exportPath,"vodaSites.txt"))
# Method #2
# the output is an sp object
vfSites <- spCbind(vfSites,joinedSitesA)
vfSites <- spCbind(vfSites,joinedSitesB)
vfSites.df <- bind_cols(vfSites@data,as.data.frame(vfSites@coords))
rm(joinedSitesA,joinedSitesB)


# Plots -------------------------------------------------------------------
# Just some not so fancy plots
tmap_mode("view")
#tmap_mode("plot")
p <- ggplot(vfSites.df,aes(x=Long,y=Lat))
p + geom_point(aes(colour=Vendor))
qplot(Long,Lat,data=vfSites.df,colour=Region)
qtm(vfSites,dots.col="Vendor",title="VF Sites")
qmplot(Long, Lat, data = vfSites.df,geom = "point") #, color = Vendor)
egBoundaries <- as.vector(spTransform(vfSites,wgs84CRS)@bbox)
names(egBoundaries) <- c("left","bottom","right","top")
egMap <- get_stamenmap(egBoundaries, maptype = "toner-lite",zoom=7)
ggmap(egMap,base_layer = p) +
  geom_point(aes(colour=Vendor),size=0.5)
tm_shape(vfSites) +
  tm_dots(col="Vendor") +
  tm_compass()


# Intersite distances metrics ---------------------------------------------

sitesCords <- vfSites@coords
k=6 # no of nearest nbrs to consider

# Metrics based on K nearest sites
distMatrix <- nn2(sitesCords,k=k)
vfSites$minDist <- apply(distMatrix[["nn.dists"]][,2:k],1,min)
vfSites$maxDist <- apply(distMatrix[["nn.dists"]][,2:k],1,max)
vfSites$meanDist <- apply(distMatrix[["nn.dists"]][,2:k],1,mean)
vfSites$sdDist <- apply(distMatrix[["nn.dists"]][,2:k],1,sd)
vfSites <- vfSites %>% mutate(normSd=sdDist/meanDist)

#temp export, for debugging
#write_csv(vfSites@data,paste0(exportPath,"vfSitestmp.csv"))

# Hierarical Clustering ---------------------------------------------------

# Hierarchical Clustering
clustFeatures <- cbind(vfSites$meanDist,vfSites$normSd)
rownames(clustFeatures) <- vfSites$Site
colnames(clustFeatures) <- c("meanDist","SD/Mean")
clustDist <- dist(clustFeatures, method = "euclidean") # distance matrix
hclust.fit <- hclust(clustDist, method="ward.D2")
# Dendogram plot
#plot(hclust.fit)
# Colored dendogram
dendSites <- as.dendrogram(hclust.fit)
dendColored <- color_branches(dendSites, k = 4)
plot(dendColored,leaflab="none") #,ylim=c(1e+05,6e+05))

# Using k=4
k=4
vfSites$SiteClass <- as.factor(cutree(hclust.fit, k=k))
levels(vfSites$SiteClass) <- c("Urban","Rural","Road","Remote")
table(vfSites$Region,vfSites$SiteClass)

#plot the clusters
qtm(vfSites,dots.col="SiteClass",title="VF Sites")

#save the sites df (in Lat/Long format)
vfSitestmp <- spTransform(vfSites,wgs84CRS)
vfSites.df <- bind_cols(vfSitestmp@data,as.data.frame(vfSitestmp@coords))
rm(vfSitestmp)
write_csv(vfSites.df,paste0(exportPath,"vodaSites_201811b.txt"))

# Further hclust trials ---------------------------------------------------

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


