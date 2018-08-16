library(readr)
library(dplyr)
library(sp)
library(tidyr)
library(rgeos)
library(rgdal)
library(spdplyr)



cellMap <- read.delim("~/R/vodaGIS/sectorMapping.txt")
cellMap$Carrier<-as.factor(cellMap$Carrier)
rownames(cellMap)<-cellMap$Suffix

atollDb <- read_delim("~/Dropbox/Voda/GISWay/sitesData/atollDb.txt","\t",
                      escape_double = FALSE, trim_ws = TRUE, skip = 1,
                      col_names = c("Site","Lat","Long"),
                      col_types = "cdd")

cellDb <- read_delim("~/Dropbox/Voda/GISWay/sitesData/oneCellDb.txt","\t",
                     escape_double = FALSE, trim_ws = TRUE,col_types = "cc")

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
    } else {
      Sector="Invalid"
    }
  return(paste(PhySite,Type,Sector,sep=","))
}

cellDb <- cellDb %>%
  rowwise %>%
  mutate(siteData=getPhySite(CellName)) %>%
  separate(siteData,into=c("Site","Type","Sector"),sep=",") %>%    
  filter(Site!="Invalid") %>%
  ungroup()

siteDb <- cellDb %>%
  group_by(Site) %>%
  summarise(noCells=n(),noSecs=n_distinct(Sector))

sites <- atollDb %>%
  inner_join(siteDb)

write_csv(sites,"~/Dropbox/Voda/GISWay/sitesData/vodaSites.txt")

coordinates(sites) <- c("Long","Lat")
proj4string(sites) <- CRS("+init=epsg:4326")
utmSites <-spTransform(sites, CRS("+init=epsg:32636"))

vodaRegions <- readOGR("/home/hnagaty/Dropbox/Voda/GISWay/geoData","Regions")
vodaRegions <-spTransform(vodaRegions, CRS("+init=epsg:32636"))
vodaRegions <- vodaRegions %>% select(REGION) %>% rename(Region=REGION)

sdg <- over(utmSites,vodaRegions)
sd <- sites@data
sdn <- bind_cols(sd,sdg)
