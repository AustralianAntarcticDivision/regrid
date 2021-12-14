## two helper functions for raster
index_coords <- function(x) {
  cell <- seq_len(raster::ncell(x))
  raster::setValues(x, cbind(raster::colFromCell(x, cell),
                             raster::rowFromCell(x, rev(cell))))
}
index_r <- function(x) {
  raster::setExtent(x, raster::extent(0, ncol(x), 0, nrow(x)))
}
## helper fun to wrap longitude into -180,180
make_180 <- function(x) {
  # lon <-   (x[,1] - 180)%%360 + 180
  # x[,1] <- lon
  lon <- x[,1]
  a360 <- lon > 180
  if (any(a360)) lon[a360] <- lon[a360] - 360
  a_360 <- lon < -180
  if (any(a_360)) lon[a_360] <- lon[a_360] + 360
  x[,1] <- lon
  x
}
## helper fun to get longlat from grid
longlatFromGrid <- function(x) {
  if (raster::isLonLat(x)) return(raster::xyFromCell(x, seq_len(raster::ncell(x))))
  xy <- xyFromCell(x, 1:raster::ncell(x))
  suppressWarnings(rgdal::project(xy, raster::projection(x), inv = TRUE))
  ## when this is fast enough
  #terra::geom(terra::project(terra::vect(xy, crs = raster::projection(x)), "OGC:CRS84"))[, c("x", "y"), drop = FALSE]
}
raster0 <- function(x, ...) {
  x <- raster::raster(x, ..., stopIfNotEqualSpaced = FALSE)
  raster::setExtent(x, raster::extent(0, ncol(x), 0, nrow(x)))
}
regrid <- function(x, coords = NULL, grid = NULL,  ...) {
  UseMethod("regrid")
}

regrid.character <- function(x, coords = NULL, grid = NULL,  ...) {
  if (!is.null(coords)) {
    if (inherits(coords, "BasicRaster")) {
      #coords <- raster::brick(raster0(x, varname = coords[1]), raster0(x, coords[2]))
    } else {
      coords <- raster::brick(raster0(x, varname = coords[1]), raster0(x, varname = coords[2]))
    }
  }
  ## here try different names
  if (is.null(coords)) {
    coords <- raster::brick(raster0(x, varname = "lon"), raster0(x, varname = "lat"))
  }
  regrid(raster0(x,...), coords = coords, grid = grid)
}
regrid.BasicRaster <- function(x, coords = NULL, grid = NULL, ...) {
  ## barycentric lookup in index-space
  if (is.null(coords)) stop("coords must be set for raster input")
  index <- quadmesh::bary_index(x, coords = index_coords(coords), grid = index_r(x))
  ok <- !is.na(index$idx)
  ## map between target grid coordinates and the index grid
  XY <- raster::values(coords)

  if (is.null(grid)) grid <- raster::raster(ncols = 360, nrows = 180)
  LL <- longlatFromGrid(grid)
  bad <- !is.finite(LL[,1]) | is.na(LL[,2])
  nn <- RANN::nn2(make_180(XY)[ok, ], LL[!bad, ], k = 1)

  iidx <- seq_len(ncell(grid))
  iidx[bad] <- NA
  r <- setValues(grid, NA_integer_)
  r[iidx[!is.na(iidx)]] <- colSums(matrix(values(x)[index$tri[, index$idx[ok]]], nrow = 3) * t(index$p)[, ok])[nn$nn.idx[,1]]
  r
}

