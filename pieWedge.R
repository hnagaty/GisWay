# function to draw a pie widge
#https://gis.stackexchange.com/questions/31294/how-to-visualize-azimuthal-data-with-uncertainties


bowtie <- function(azimuth, delta, origin=c(0,0), radius=1, eps=1) {
  #
  # On entry:
  #   azimuth and delta are in degrees.
  #   azimuth is east of north; delta should be positive.
  #   origin is (lon, lat) in degrees.
  #   radius is in meters.
  #   eps is in degrees: it is the angular spacing between vertices.
  #
  # On exit:
  #   returns an n by 2 array of (lon, lat) coordinates describing a "bowtie" shape.
  #
  # NB: we work in radians throughout, making conversions from and to degrees at the
  #   entry and exit.
  #--------------------------------------------------------------------------------#
  if (eps <= 0) stop("eps must be positive")
  if (delta <= 0) stop ("delta must be positive")
  if (delta > 90) stop ("delta must be between 0 and 90")
  if (delta >= eps * 10^4) stop("eps is too small compared to delta")
  if (origin[2] > 90 || origin[2] < -90) stop("origin must be in lon-lat")
  a <- azimuth * pi/180; da <- delta * pi/180; de <- eps * pi/180 
  start <- origin * pi/180
  #
  # Precompute values for `goto`.
  #
  lon <- start[1]; lat <- start[2]
  lat.c <- cos(lat); lat.s <- sin(lat)
  radius.radians <- radius/6366710
  radius.c <- cos(radius.radians); radius.s <- sin(radius.radians) * lat.c
  #
  # Find the point at a distance of `radius` from the origin at a bearing of theta.
  # http://williams.best.vwh.net/avform.htm#Math
  #
  goto <- function(theta) {
    lat1 <- asin(lat1.s <- lat.s * radius.c + radius.s * cos(theta))
    dlon <- atan2(-sin(theta) * radius.s, radius.c - lat.s * lat1.s)
    lon1 <- lon - dlon + pi %% (2*pi) - pi
    c(lon1, lat1)
  }
  #
  # Compute the perimeter vertices.
  #
  n.vertices <- ceiling(2*da/de)
  bearings <- seq(from=a-da, to=a+da, length.out=n.vertices)
  t(cbind(start,
          sapply(bearings, goto),
          start,
          sapply(rev(bearings+pi), goto),
          start) * 180/pi)
}

n <- 10
input <- data.frame(cbind(
  id = 1:n, 
  lon = runif(n, -180, 180),
  lat = asin(runif(n)) * 180/pi,
  azimuth = runif(n, 0, 360),
  delta = 90 * rbeta(n, 20, 70),
  radius = 10^7/90 * rgamma(n, 10, scale=2/10)
))

shapes <- as.data.frame(do.call(rbind, 
                                by(input, input$id, 
                                   function(d) cbind(d$id, bowtie(d$azimuth, d$delta, c(d$lon, d$lat), d$radius, 1)))))
colnames(shapes) <- c("id", "x", "y")
plot(shapes$x, shapes$y, type="n", xlab="Longitude", ylab="Latitude", main="Bowties")
temp <- by(shapes, shapes$id, function(d) lines(d$x, d$y, type="l", lwd=2, col=d$id))
