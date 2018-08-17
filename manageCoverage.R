library(readr)
library(sp)
library(rgdal)
library(raster)
library(rgeos)

cellCoverage <- readOGR("D:/Optimisation/~InProgress/201806_GisFramework/pilotCoverage.geojson", "pilotCoverage")
siradel <- raster("D:/Optimisation/~InProgress/201806_GisFramework/ClippedDLU")
siradelLegend <- read_delim("D:/Optimisation/~InProgress/201806_GisFramework/geoData/Siradel2016_DLU/EGYPT_3_DLU_50m.mnu",
                      " ",col_names = c("Code","Clutter"))


siradel <- as.factor(siradel)
plot(siradel)


intersectingCoverage <- gIntersection(cellCoverage,cellCoverage,byid=TRUE)

test <- gUnaryUnion(cellCoverage)

test2 <- rasterize(test,siradel,fun=sum)
