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

lasaggregatepolygons = function(las, type, concavity, length_threshold, by)
{
  X <- Y <- tree <- NULL

  if (type == "convex")
    dt = las@data[, if (!anyNA(.BY)) stdtreehullconvex(X,Y, .GRP), by = by]
  else
    dt = las@data[, if (!anyNA(.BY)) stdtreehullconcave(X,Y, .GRP, concavity, length_threshold), by = by]

  data.table::setnames(dt, names(dt), c(by, "polygon"))

  polygons = dt$polygon
  polygons = polygons[!sapply(polygons, is.null)]

  if (length(polygons) == 0)
    return(NULL)

  spoly = sp::SpatialPolygons(polygons)

  for (i in 1:length(spoly)) spoly@polygons[[i]]@ID = as.character(i)

  data = data.frame(dt[, 1])
  spdf = sp::SpatialPolygonsDataFrame(spoly, data)
  sp::proj4string(spdf)<-las@crs

  return(spdf)
}

stdtreehullconvex = function(x,y, grp)
{
  if (length(x) < 4)
    return(NULL)

  i = grDevices::chull(x,y)
  i = c(i, i[1])
  P = cbind(x[i], y[i])
  poly = sp::Polygon(P)
  poly = sp::Polygons(list(poly), ID = grp)

  list(poly = list(poly))
}

stdtreehullconcave = function(x,y, grp, concavity, length_threshold)
{
  if (length(x) < 4)
    return(NULL)

  P = concaveman::concaveman(cbind(x,y), concavity, length_threshold)
  poly = sp::Polygon(P)
  poly = sp::Polygons(list(poly), ID = grp)

  list(poly = list(poly))
}

