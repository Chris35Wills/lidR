# ===============================================================================
#
# PROGRAMMERS:
#
# jean-romain.roussel.1@ulaval.ca  -  https://github.com/Jean-Romain/rlas
#
# COPYRIGHT:
#
# Copyright 2018 Jean-Romain Roussel
#
# This file is part of lidR R package.
#
# rlas is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#
# ===============================================================================

#' Find the flighlines hulls
#'
#' Find the flightline hulls from a catalog. The function works only if the field 'PointSourceID'
#' is properly populated. The function first loads a random fraction of the catalog then compute a concave
#' hull for each flightline. Flighlines are expected to be recorded in 'PointSourceID'. If not, the function
#' will run anyway but will find a single flighline being the whole dataset itself.
#'
#' The flightlines should be contiguous. You could get weird results when the flightline IDs are spread
#' in several non contiguous parts of the point cloud. In other words "islands" of points are not
#' supported yet.\cr\cr
#' The concave hull method under the hood is described in Park & Oh (2012). The function relies on
#' the \link[concaveman:concaveman]{concaveman} function which itself is a wrapper around the
#' \href{https://github.com/mapbox/concaveman}{Vladimir Agafonking's implementation}.\cr\cr
#'
#' @param ctg A \link[lidR:catalog]{LAScatalog} object.
#' @param concavity numeric. If \code{type = "concave"}, a relative measure of concavity. 1 results
#' in a relatively detailed shape, Infinity results in a convex hull.
#' @param length_threshold numeric. If \code{type = "concave"}, when a segment length is under this
#' threshold, it stops being considered for further detalization. Higher values result in simpler shapes.
#' @param density numeric. To compute the concave hull of each flighline in a reasonnable timeframe and using
#' a little amount of memory, only a random fraction of the point is loaded. This number represents the
#' number of points point per square unit actually loaded. Default is 0.2 i.e. approximately 0.2 pt per
#' square meter. This is enought to get a fairly accurate hulls (~ 1 point each 2 meters).
#' @return A \code{SpatialPolygonDataFrame}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using files for which the PointSourceID is properly populated
#' ctg = catalog("folder/")
#' by_file(ctg) = TRUE
#' buffer(ctg) = 5
#'
#' flightlines =  catalog_flightlines(ctg)
#'
#' col = adjustcolor(col = pastel.colors(length(flightlines)), alpha.f = 0.6)
#' sp::plot(flightlines, col = col)
#' }
catalog_flightlines = function(ctg, concavity = 2, length_threshold = 50, density = 0.2)
{
  message("This is un experimental function. User should check twice the output!")

  npoints  <- sum(ctg@data$`Number of point records`)

  progress <- progress(ctg)
  ncores   <- cores(ctg)
  surface  <- area(ctg)

  pdensity <- npoints/surface
  fraction <- round(density/pdensity, 3)

  filter   <- paste("-keep_random_fraction", fraction)

  clusters <- catalog_makecluster(ctg, 1)
  nclust   <- length(clusters)

  if (nclust < ncores)
    ncores <- nclust

  future::plan(future::multiprocess, workers = ncores)

  spdf = list()
  for(i in seq_along(clusters))
  {
    cluster = clusters[[i]]

    spdf[[i]] <- future::future(
    {
      las = suppressMessages(readLAS(cluster, select = "xyzp", filter = filter))

      if (is.null(las))
        return(NULL)

      lasaggregatepolygons(las, "concave", concavity, length_threshold, by = "PointSourceID")
    }, earlySignal = TRUE)

    if(progress)
    {
      cat(sprintf("\rProgress: %g%%", round(i/nclust*100)), file = stderr())
      graphics::rect(cluster@bbox$xmin, cluster@bbox$ymin, cluster@bbox$xmax, cluster@bbox$ymax, border = "black", col = "forestgreen")
    }
  }

  if (progress)
    message("\nCatalog processed. Dissolving the polygons... This may take several seconds.")

  spdf = future::values(spdf)
  spdf = spdf[!sapply(spdf, is.null)]
  spdf = do.call(sp::rbind.SpatialPolygonsDataFrame, spdf)
  spdf = rgeos::gUnaryUnion(spdf, spdf@data[,1])

  return(spdf)
}
