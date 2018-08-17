library("sp")
library(readxl)
library(SearchTrees)
library(rgeos)

Sites <- read_excel("D:/Optimisation/~InProgress/201806_GisFramework/grbSites.xlsx")
coordinates(Sites) <- c("Longitude", "Latitude")
summary(Sites)
crs <- CRS("+init=EPSG:32636")
proj4string(Sites) <- crs
summary(Sites)

SitesC <- Sites
tree <- createTree(coordinates(Sites))
inds <- knnLookup(tree, newdat = coordinates(SitesC), k = 6)

avgDist <- function(sites) {
  return(mean(gDistance(sites,byid=TRUE)[2:6]))
}

knnsites <- Sites[inds[1,],1]
avgDist(knnsites)


