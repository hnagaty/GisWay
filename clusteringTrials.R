# Ward Hierarchical Clustering
d <- dist(sitesCl, method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
#rect.hclust(fit, k=5, border="red")
sitesDf$Cluster <- as.factor(groups)
ggplot(sitesDf,aes(x=Long,y=Lat,col=Cluster)) + geom_point()

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