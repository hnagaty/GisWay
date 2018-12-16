
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(sp)

# Shared variables --------------------------------------------------------

#myDesktop <- "Linux"
myDesktop <- "Windows"

if (myDesktop == "Linux") {
  paths <- list()
  paths$desktop <- myDesktop
  paths$parentPath <- "~/Dropbox/Voda/GISWay/"
  paths$dataDir <- "~/data/gmrr/"
  paths$exportPath <- "~/Dropbox/Voda/GISWay/export/"
  paths$pngPath <- "~/Dropbox/Voda/GISWay/gisWay/outs/plots/"
  paths$localPath <- "~/data/geodata/Siradel_2016March/DLU/"
  paths$outDir <- "~/Dropbox/Voda/GISWay/export/"
} else {
  paths <- list()
  paths$desktop <- myDesktop
  paths$parentPath <- "D:/Optimisation/~InProgress/201806_GisFramework/"
  paths$dataDir<- "d:/data/mrr/2018Sep25/"
  paths$exportPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
  paths$pngPath <- "D:/Optimisation/~InProgress/201806_GisFramework/gisWay/outputs/plots/"
  paths$geoPath <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
  paths$ossDataDir <- "D:/Optimisation/~InProgress/201806_GisFramework/ossData/"
  paths$miscDataDir <- "D:/Optimisation/~InProgress/201806_GisFramework/miscData/"
  paths$localPath <- "D:/Optimisation/~InProgress/201806_GisFramework/geoData/Siradel_2016_Country/DLU/"
  paths$outDir <- "D:/Optimisation/~InProgress/201806_GisFramework/export/"
  paths$wmrrDir <- "D:/Optimisation/~InProgress/201806_GisFramework/ossData/wmrr_122018_dec03/"
  paths$gmrrDir <- "D:/Optimisation/~InProgress/201806_GisFramework/ossData/mrr_092018_sep25/"
}

rm(myDesktop)

paths$myPath <-  paths$dataDir
paths$mainPath <- paste0(paths$parentPath,"gisWay/")
paths$globalPath <- paste0(paths$parentPath,"sitesData/")
paths$exportPath <- paste0(paths$parentPath,"export/")
paths$geoPath <- paste0(paths$parentPath,"geoData")

files <- list()
files$atollFile <- "atollSites.csv"
files$cellFile <- "oneCell.csv"

egUtmCRS <- CRS("+init=epsg:32636")
wgs84CRS <- CRS("+init=epsg:4326")


# Shared Functions --------------------------------------------------------

# Needed for GMRR parsing
pasteDir <- function(c) {
  if (nchar(c) == 2) {return (c)}
  else if (c=="BSPOWER") {return("DL")}
  else if (c=="MSPOWER") {return("UL")}
  else {return("BL")}
}

# Scans a folder and checks for complete MRR files in that folder
getMrrFiles <- function(dataDir) {
  allMrrFiles <- list.files(dataDir,pattern=".conf")
  allMrrFiles <- sub(".conf","",allMrrFiles)
  
  mrrConf.df <- data.frame(File=as.character(),
                           Creator=as.character(),
                           StartTime=as.character(),
                           EndTime=as.character(),
                           noCells=as.numeric(),
                           Exists=as.logical())
  
  for (mrrFile in allMrrFiles) {
    mrrConf <- read_lines(paste0(dataDir,mrrFile,".conf"))
    mrrStartTime <-  mrrConf[grepl("start time", tolower(mrrConf))]
    mrrStartTime <- substring(mrrStartTime,regexpr("\t",mrrStartTime)+1)
    mrrEndTime <-  mrrConf[grepl("stop time", tolower(mrrConf))]
    mrrEndTime <- substring(mrrEndTime,regexpr("\t",mrrEndTime)+1)
    mrrCreator <-  mrrConf[grepl("creator", tolower(mrrConf))]
    mrrCreator <- substring(mrrCreator,regexpr("\t",mrrCreator)+1)
    mrrNoCells <-  mrrConf[grepl("no of cells", tolower(mrrConf))]
    mrrNoCells <- as.numeric(substring(mrrNoCells,regexpr("\t",mrrNoCells)+1))
    mrrExists <- file.exists(paste0(dataDir,mrrFile,".msmt"))
    mrr.df <- data.frame(File=mrrFile,
                         Creator=mrrCreator,
                         StartTime=mrrStartTime,
                         EndTime=mrrEndTime,
                         noCells=mrrNoCells,
                         Exists=mrrExists)
    mrrConf.df <- rbind.data.frame(mrrConf.df,mrr.df)
  }
 
  mrrConf.df$StartTime <- as.POSIXct(mrrConf.df$StartTime)
  mrrConf.df$EndTime <- as.POSIXct(mrrConf.df$EndTime)
  mrrConf.df <- mrrConf.df %>%
    filter(Exists=="TRUE") %>%
    arrange(desc(StartTime))
  return(mrrConf.df)
}


# Read mrr files in a given folder. Reads the files mentioned in the df outputed from getMrrFiles()
readMrrFiles <- function(path, filesDf) {
  mrrFiles <- paste0(filesDf$File,".msmt")
  df <- data.frame()
  for (m in mrrFiles) {
    mrrtmp <- read_tsv(paste0(path,m))
    df <- bind_rows(df,mrrtmp)
  }
  return(df)
}

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


# hClustMrr function
hClustMrr <- function (data, labels, infoList, k=2, save=TRUE, plotDendogram=TRUE) {
  scaled <- scale(data)
  fit <- hclust.vector(scaled,method="ward",metric="euclidean")
  fit$info <- cList
  # Dendogram plot
  if (plotDendogram) {
    print("Now plotting the dendogram ...")
    dendMrr <- as.dendrogram(fit)
    dendColored <- color_branches(dendMrr, k = k) #,groupLabels = TRUE) ==> seems to not match clusters in cuttree
    plot(dendColored,leaflab="none")
    if (save) {
      dev.print(png, paste0(paths$exportPath,infoList$class, "_HClustDendrogram_",infoList$version,".png"), width=1600, height=900)
    }
  }
  # Cut the tree
  clusters <- as.factor(cutree(fit, k = k))
  print(table(clusters))
  clustersTable <-  bind_cols(labels,data.frame(cluster=clusters))
  if (save) {
    clstList <- list(info=infoList,clusters=clustersTable)
    write_lines(clstList$info,path = paste0(paths$exportPath, infoList$class, "_HClustOutput_",infoList$version,"_info.txt"), sep="\r\n")
    write_lines(table(clusters),path = paste0(paths$exportPath, infoList$class, "_HClustOutput_",infoList$version,"_info.txt"), sep=" ", append = TRUE)
    write_csv(clstList$clusters, path = paste0(paths$exportPath,infoList$class, "_HClustOutput_",infoList$version,"_clusters.csv"))
    saveRDS(clstList, file = paste0(paths$exportPath, infoList$class, "_HClustOutput_",infoList$version,".rds"))
  }
  returnList <- (list(model=fit,clusters=clustersTable))
  if (save) {
    modelFilename <- paste0(paths$exportPath,infoList$class, "_HClustFit_",infoList$version,".rds")
    saveRDS(returnList,file=modelFilename)
    message(paste0("Model is saved in ", modelFilename))
  }
  return(returnList)
}

tryHClust <- function(model, labeled, h, k) {
  clstr <- cutree(model, k = k, h = h)
  print(table(clstr))
  print("Plotting the dendrogram ...")
  dend <- as.dendrogram(model)
  #order <- order.dendrogram(dend)
  dendColored <- color_branches(dend, h = h, k = k) #groupLabels = TRUE,  clusters = clstr[order]) #But the labels are not printed
  plot(dendColored,leaflab="none")
  title(main=model$info$version, sub=paste(sep = "\n", model$info$notes, paste0("h=",h,", k=",k)))
  labeled$cluster <- as.factor(clstr)
  return(labeled)
}

#remove spaces and other characters from kpi names
strpKpiName <- function(x) { 
  tolower(gsub("[ /]","",x))  
} 

# Read shared data tables -------------------------------------------------

# The WMRR bin mapping
binMap <- read_tsv("wmrrRanges.txt",
                   skip = 1,
                   col_names = c("Bin", "Bler", "UeTxPwr", "CpichEcno", "CpichRscp","DlTxCodePowerSf"),
                   col_types = "iccccc")
# Cellmap file, maps cell suffix to sector & carrier. This is not the latest version
cellMap <- read.delim(paste0(paths$mainPath,"sectorMapping.txt"))
cellMap$Carrier<-as.factor(cellMap$Carrier)
rownames(cellMap)<-cellMap$Suffix
