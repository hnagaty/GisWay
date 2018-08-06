library(readr)
library(dplyr)


cellMap <- read_delim("sectorMapping.txt", "\t",
    escape_double = FALSE, na = "NA",
    col_types = cols(Carrier = col_factor(levels = c("1", "2", "3")),
                     Sector = col_factor(levels = c("S1", "S2", "S3", "S4", "S5"))),
    trim_ws = TRUE)

cellMap <- read.delim("~/R/vodaGIS/sectorMapping.txt")
cellMap$Carrier<-as.factor(cellMap$Carrier)
cellMap$rownames<-cellMap$Suffix


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
    {
      a
    } else {
      Sector="Invalid"
    }
  }
  return(c(PhySite,Type,Sector))
}

cellDb <- cellDb %>%
  rowwise %>%
  mutate(Site=getPhySite(CellName)[1],Type=getPhySite(CellName)[2]) %>%
  filter(Site!="Invalid") %>%
  ungroup()

siteDb <- cellDb %>%
  group_by(Site) %>%
  summarise(noCells=n())

allDb <- atollDb %>%
  inner_join(siteDb)
