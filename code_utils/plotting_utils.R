##### DENSITY COLORED SCATTER PLOTS #####

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Define a stat_scatter_density helper for ggplot2
StatDensity <- ggproto("StatDensity", Stat,
  compute_group = function(data, scales) {
    data$density <- get_density(data$x, data$y)
    data
  },
  
  required_aes = c("x", "y"),
  default_aes = aes(size=1)
)

#' Scatter data and color by density. Works like geom_point, but adds a density
#' variable and colors the points by density
stat_scatter_density <- function(mapping = aes(), data = NULL, geom = "point",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  # A little magic to overwrite the color parameter by default
  mapping$colour <-  quote(stat(density))
  layer(
    stat = StatDensity, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
