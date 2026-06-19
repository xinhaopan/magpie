# Extract N Load (nitrogen surplus) from MAgPIE GDX for 2020 and 2050
# and write NetCDF files into the corresponding SEALS output directories.

library(magpie4)
library(magclass)
library(raster)
library(ncdf4)

# writeRasterBrick helper (adapted from disaggregation_LUH2.R)
writeRasterBrick <- function(x, filePath, comment = NULL, zname = "Time", ...) {
  tmp <- names(x)
  tmp <- strsplit(tmp, "\\.\\.")
  years    <- sort(unique(unlist(lapply(tmp, function(x) x[1]))))
  varnames <- sort(unique(unlist(lapply(tmp, function(x) x[2]))))
  zunit <- ifelse(all(isYear(years)), "years", "")
  years <- as.numeric(gsub("y", "", years))
  if (is.null(varnames)) varnames <- "Variable"
  if (is.null(comment)) comment <- "not specified"

  nc_dims <- list(
    ncdf4::ncdim_def("lon", "degrees_east",  sort(unique(raster::xFromCell(x, seq_len(raster::ncell(x)))))),
    ncdf4::ncdim_def("lat", "degrees_north", sort(unique(raster::yFromCell(x, seq_len(raster::ncell(x)))))),
    ncdf4::ncdim_def(zname, zunit, years)
  )
  nc_vars <- lapply(varnames, function(v) ncdf4::ncvar_def(v, comment, nc_dims, -9999))
  nc_file <- ncdf4::nc_create(filePath, nc_vars)
  for (v in varnames) {
    arr <- raster::as.array(x[[grep(v, names(x))]])
    arr <- aperm(arr, c(2, 1, 3))
    arr[is.na(arr)] <- -9999
    ncdf4::ncvar_put(nc_file, v, arr)
  }
  ncdf4::nc_close(nc_file)
  message("Written: ", filePath)
}

base <- "/home/xinhaopa/work/projects/MAgPIE-SEAL/magpie/output"

scenarios <- list(
  ref = list(
    gdx_dir  = file.path(base, "Global_Ref2050"),
    out_dirs = c(
      file.path(base, "seals_europe_ref"),
      file.path(base, "seals_indonesia_ref")
    )
  ),
  ghg = list(
    gdx_dir  = file.path(base, "Global_GHG2050"),
    out_dirs = c(
      file.path(base, "seals_europe_ghg"),
      file.path(base, "seals_indonesia_ghg")
    )
  )
)

years <- c("y2020", "y2050")

for (scen_name in names(scenarios)) {
  scen    <- scenarios[[scen_name]]
  gdx_dir <- scen$gdx_dir
  gdx     <- file.path(gdx_dir, "fulldata.gdx")

  message("\n=== Processing: ", scen_name, " ===")
  nb      <- NitrogenBudget(gdx, level = "grid")
  surplus <- nb[, years, "surplus"]
  surplus <- clean_magpie(collapseNames(surplus, collapsedim = 3.1))

  # Attach spatial dimension metadata so as.RasterBrick works
  getSets(surplus, fulldim = FALSE)[1] <- "x.y.iso"

  message("Converting to raster...")
  r <- as.RasterBrick(surplus)

  for (out_dir in scen$out_dirs) {
    nc_path <- file.path(out_dir, paste0("NLoad_", scen_name, "_2020_2050.nc"))
    message("Writing: ", nc_path)
    raster::writeRaster(r, nc_path, format = "CDF", overwrite = TRUE,
                        varname = "NLoad", varunit = "kgN_per_ha",
                        longname = "Nitrogen surplus (N Load)",
                        xname = "lon", yname = "lat", zname = "time")
    message("Done.")
  }
}

message("\nAll N Load files written successfully.")
