# Supplementary code to calculate the average intersite distance
# It's now incorporated into readSites.R
# This script also has the clustering trials

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



egUtmCRS <- CRS("+init=epsg:32636")
wgs84CRS <- CRS("+init=epsg:4326")

sites <- read_delim("D:/Optimisation/~InProgress/201806_GisFramework/gisWay/atollSampleSites.txt","\t",
                      escape_double = FALSE, trim_ws = TRUE, skip = 1,
                      col_names = c("Site","Lat","Long"),
                      col_types = "cdd")
# Remove duplicate sites
sites = sites[order(sites$Site),]
sites = sites[!duplicated(sites$Site),]
#write_csv(sites,"vodaSites.txt")

coordinates(sites) <- c("Long","Lat")
proj4string(sites) <- wgs84CRS
utmSites <-spTransform(sites, egUtmCRS)
sitesDf <- bind_cols(utmSites@data,as.data.frame(utmSites@coords))

sitesCord <- sitesDf %>%
  dplyr::select(Long,Lat)

sitesNames <- sitesDf$Site

k=6
noSites <- NROW(sitesNames)

# not used library
#sitesDistances <- distances(sitesDf,id_variable = "Site",dist_variables = c("Long","Lat"))
#nearestSitesIdx <- nearest_neighbor_search(sitesDistances,k=k)
#nearestSitesNames <- matrix(sitesNames[nearestSitesIdx],nrow=k,ncol=noSites)

distMatrix <- nn2(sitesCord,k=k)
sitesDf$minDist <- apply(distMatrix[["nn.dists"]][,2:k],1,min)
sitesDf$maxDist <- apply(distMatrix[["nn.dists"]][,2:k],1,max)
sitesDf$meanDist <- apply(distMatrix[["nn.dists"]][,2:k],1,mean)
sitesDf$sdDist <- apply(distMatrix[["nn.dists"]][,2:k],1,sd)
sitesDf <- sitesDf %>%
  mutate(normSd=sdDist/meanDist)

#write_csv(sites,"vodaSitesDf.txt")

# Clustering
sitesCl <- sitesDf %>%
  dplyr::select(meanDist,normSd)

# Hierarchical Clustering
hclust.dist <- dist(sitesCl, method = "euclidean") # distance matrix
hclust.fit <- hclust(hclust.dist, method="ward.D2")
# Single k value
#plot(hclust.fit) # display dendogram
#sites.clusters <- cutree(hclust.fit, k=2) # cut tree into k clusters
#sitesDf$Cluster <- as.factor(sites.clusters)
#p <- ggplot(data=sitesDf,aes(x=Long,y=Lat)) + geom_point()
#p + aes(col=Cluster)
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

p1 <- ggplot(sitesCl,aes(x=meanDist)) + geom_histogram()
p1

sites.export <- dplyr::filter(sites.clusters,kValue==5)
write_csv(sites.export,"clusteredSites.txt")
# draw dendogram with red borders around the 5 clusters
#rect.hclust(fit, k=2, border="red") 


# kmeans
# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
# K-Means Cluster Analysis
fit <- kmeans(mydata, 2) # 5 cluster solution
# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)
