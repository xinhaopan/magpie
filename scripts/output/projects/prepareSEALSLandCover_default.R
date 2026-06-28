# Convert regional SEALS land-cover outputs back to ESA CCI-style land-cover
# codes and collect the matching NLoad NetCDF files.

library(terra)
library(ncdf4)
library(foreach)
library(doParallel)

projectName <- "GlobalLandUse"
sealsBaseYr <- 2020

# The script is meant to be run from either the repository root or magpie/.
repoDir <- normalizePath(getwd(), mustWork = TRUE)
magpieDir <- if (basename(repoDir) == "magpie") repoDir else file.path(repoDir, "magpie")
if (!dir.exists(magpieDir)) {
  stop("Could not find the magpie directory. Run this script from the repo root or magpie/.")
}

dirMagOut <- file.path(magpieDir, "output")
dirOut <- file.path(dirMagOut, "seals_landcover_esacci")
dirTmp <- file.path(tempdir(), paste0("prepareSEALSLandCover_", projectName))

# ESA CCI reference file with original land-cover codes. Unchanged SEALS pixels
# keep these original ESA codes; changed pixels are mapped from SEALS7 to ESA-like
# codes below.
esaBaseFile <- "/p/projects/magpie/users/vjeetze/data/ESA_CCI_land_cover/product/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.tif"

regions <- list(
  europe = list(
    baseline = file.path(dirMagOut, "seals_baseline_2020", "europe",
                         "lulc_esa_seals7_Baseline_2020_clipped.tif"),
    scenarios = list(
      ref = file.path(dirMagOut, "seals_europe_ref"),
      ghg = file.path(dirMagOut, "seals_europe_ghg")
    )
  ),
  indonesia = list(
    baseline = file.path(dirMagOut, "seals_baseline_2020", "indonesia",
                         "lulc_esa_seals7_Baseline_2020_clipped.tif"),
    scenarios = list(
      ref = file.path(dirMagOut, "seals_indonesia_ref"),
      ghg = file.path(dirMagOut, "seals_indonesia_ghg")
    )
  )
)

seals7ToEsacci <- rbind(
  c(1, 190), # urban
  c(2, 10),  # cropland
  c(3, 130), # grassland
  c(4, 60),  # forest
  c(5, 125), # other natural / non-forest vegetation
  c(6, 210), # water
  c(7, 200)  # bare
)

dir.create(dirOut, recursive = TRUE, showWarnings = FALSE)
dir.create(dirTmp, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = dirTmp, todisk = TRUE, memmax = 16, memmin = 1e-9)

if (!file.exists(esaBaseFile)) {
  stop("ESA CCI reference file not found: ", esaBaseFile)
}

findSealsProjectionFiles <- function(sealsDir) {
  dirIn <- file.path(sealsDir, "intermediate", "stitched_lulc_simplified_scenarios")
  if (!dir.exists(dirIn)) {
    stop("SEALS stitched scenario directory not found: ", dirIn)
  }

  files <- list.files(dirIn, pattern = "\\.tif$", full.names = TRUE)
  files <- files[!grepl(paste0(sealsBaseYr, "(_clipped)?\\.tif$"), basename(files))]

  clipped <- files[grepl("_clipped\\.tif$", basename(files))]
  if (length(clipped) > 0) {
    files <- clipped
  }

  if (length(files) == 0) {
    stop("No projected SEALS GeoTIFF files found in: ", dirIn)
  }
  files
}

makeEsacciBaseline <- function(regionName, baselineSeals7File) {
  regionOut <- file.path(dirOut, regionName)
  dir.create(regionOut, recursive = TRUE, showWarnings = FALSE)

  lcBaseFile <- file.path(regionOut, paste0("lulc_Baseline_", sealsBaseYr, "_esacci_code.tif"))
  if (file.exists(lcBaseFile)) {
    return(lcBaseFile)
  }

  message("Creating ESA CCI baseline for ", regionName)
  esaBase <- rast(esaBaseFile)
  sealsBase <- rast(baselineSeals7File)
  lcBase <- crop(esaBase, sealsBase, snap = "near")

  if (!compareGeom(lcBase, sealsBase, stopOnError = FALSE)) {
    lcBase <- resample(lcBase, sealsBase, method = "near")
  }

  writeRaster(lcBase, lcBaseFile,
              overwrite = TRUE, filetype = "GTiff",
              datatype = "INT2U", gdal = c("COMPRESS=DEFLATE"))
  lcBaseFile
}

buildTasks <- function() {
  tasks <- list()
  for (regionName in names(regions)) {
    region <- regions[[regionName]]
    if (!file.exists(region$baseline)) {
      stop("Regional SEALS7 baseline not found: ", region$baseline)
    }

    lcBaseFile <- makeEsacciBaseline(regionName, region$baseline)
    for (scenarioName in names(region$scenarios)) {
      sealsDir <- region$scenarios[[scenarioName]]
      for (projectionFile in findSealsProjectionFiles(sealsDir)) {
        outFile <- file.path(
          dirOut, regionName,
          sub("\\.tif$", "_esacci_code.tif", basename(projectionFile))
        )
        tasks[[length(tasks) + 1]] <- list(
          region = regionName,
          scenario = scenarioName,
          sealsDir = sealsDir,
          baselineSeals7 = region$baseline,
          baselineEsacci = lcBaseFile,
          projection = projectionFile,
          output = outFile
        )
      }
    }
  }
  tasks
}

copyNLoadFiles <- function() {
  copied <- character()
  for (regionName in names(regions)) {
    regionOut <- file.path(dirOut, regionName)
    dir.create(regionOut, recursive = TRUE, showWarnings = FALSE)

    for (scenarioName in names(regions[[regionName]]$scenarios)) {
      src <- file.path(regions[[regionName]]$scenarios[[scenarioName]],
                       paste0("NLoad_", scenarioName, "_2020_2050.nc"))
      if (!file.exists(src)) {
        warning("NLoad file not found: ", src)
        next
      }
      nc <- ncdf4::nc_open(src)
      nloadUnit <- ncdf4::ncatt_get(nc, "NLoad", "units")$value
      ncdf4::nc_close(nc)
      if (!identical(nloadUnit, "kgN_per_ha_per_yr")) {
        warning("Unexpected NLoad unit in ", src, ": ", nloadUnit,
                ". Regenerate it with export_NLoad_to_SEALS.R before using it.")
      }
      dst <- file.path(regionOut, basename(src))
      file.copy(src, dst, overwrite = TRUE)
      copied <- c(copied, dst)
    }
  }
  copied
}

writeReadme <- function(tasks, nloadFiles) {
  readme <- file.path(dirOut, "README.md")
  lines <- c(
    "# SEALS Land Cover And NLoad Outputs",
    "# SEALS 土地利用和 NLoad 输出",
    "",
    "This folder is generated by `scripts/output/projects/prepareSEALSLandCover_default.R`.",
    "本文件夹由 `scripts/output/projects/prepareSEALSLandCover_default.R` 生成。",
    "",
    "## Purpose / 目的",
    "",
    "The script converts regional SEALS stitched land-cover outputs from the simplified SEALS7 classes back to ESA CCI-style land-cover codes.",
    "这个脚本把区域 SEALS stitched 土地利用输出从简化的 SEALS7 类别转换回 ESA CCI-style 土地利用代码。",
    "Pixels without land-cover change keep the original ESA CCI 2020 code. Pixels changed by SEALS are reclassified from SEALS7 to ESA-style codes.",
    "没有土地利用变化的像元保留原始 ESA CCI 2020 代码；被 SEALS 改变的像元则从 SEALS7 重新分类为 ESA-style 代码。",
    "",
    "## Regions And Scenarios / 区域和情景",
    "",
    "- Regions: `europe`, `indonesia`.",
    "- 区域: `europe`, `indonesia`。",
    "- Scenarios: `ref` (`Global_Ref2050`) and `ghg` (`Global_GHG2050`).",
    "- 情景: `ref` (`Global_Ref2050`) 和 `ghg` (`Global_GHG2050`)。",
    "- Baseline year: `2020`.",
    "- 基准年: `2020`。",
    "",
    "## SEALS7 To ESA-Style Mapping / SEALS7 到 ESA-style 的映射",
    "",
    "| SEALS7 | ESA-style code | Meaning / 含义 |",
    "|---:|---:|---|",
    "| 1 | 190 | urban / 城市 |",
    "| 2 | 10 | cropland / 耕地 |",
    "| 3 | 130 | grassland / 草地 |",
    "| 4 | 60 | forest / 森林 |",
    "| 5 | 125 | other natural / non-forest vegetation / 其它自然或非森林植被 |",
    "| 6 | 210 | water / 水体 |",
    "| 7 | 200 | bare / 裸地 |",
    "",
    "## NLoad Files / NLoad 文件",
    "",
    "The matching NLoad NetCDF files are copied into each regional output folder.",
    "匹配的 NLoad NetCDF 文件会被复制到每个区域输出文件夹中。",
    "Their unit is `kgN_per_ha_per_yr`. In `export_NLoad_to_SEALS.R`, MAgPIE absolute nitrogen surplus (`mio. t N/yr`) is divided by grid-level cropland area (`mio. ha`) and multiplied by 1000.",
    "单位是 `kgN_per_ha_per_yr`。在 `export_NLoad_to_SEALS.R` 中，MAgPIE 的绝对氮盈余 (`mio. t N/yr`) 会除以 grid-level cropland area (`mio. ha`) 并乘以 1000。",
    "",
    "## Generated Files / 生成的土地利用文件",
    ""
  )

  taskLines <- vapply(tasks, function(task) {
    paste0("- `", normalizePath(task$output, mustWork = FALSE), "`")
  }, character(1))
  nloadLines <- vapply(nloadFiles, function(path) {
    paste0("- `", normalizePath(path, mustWork = FALSE), "`")
  }, character(1))

  writeLines(c(lines, taskLines, "", "## Copied NLoad Files / 复制的 NLoad 文件", "", nloadLines), readme)
  readme
}

processTask <- function(task) {
  if (file.exists(task$output)) {
    message(task$output, " already exists.")
    return(task$output)
  }

  taskTmp <- file.path(dirTmp, paste0(task$region, "_", task$scenario, "_", Sys.getpid()))
  dir.create(taskTmp, recursive = TRUE, showWarnings = FALSE)
  terraOptions(tempdir = taskTmp, todisk = TRUE, memmax = 8, memmin = 1e-9)

  lcBase <- rast(task$baselineEsacci)
  sealsBase <- rast(task$baselineSeals7)
  sealsProj <- rast(task$projection)

  if (!compareGeom(sealsBase, sealsProj, stopOnError = FALSE)) {
    sealsProj <- crop(sealsProj, sealsBase, snap = "near")
    if (!compareGeom(sealsBase, sealsProj, stopOnError = FALSE)) {
      sealsProj <- resample(sealsProj, sealsBase, method = "near")
    }
  }

  if (!compareGeom(lcBase, sealsProj, stopOnError = FALSE)) {
    lcBase <- crop(lcBase, sealsProj, snap = "near")
    if (!compareGeom(lcBase, sealsProj, stopOnError = FALSE)) {
      lcBase <- resample(lcBase, sealsProj, method = "near")
    }
  }

  lcDiff <- classify(sealsBase - sealsProj,
                     rbind(c(-Inf, 0, 1), c(1, Inf, 1)),
                     right = FALSE)
  sealsDiff <- terra::mask(sealsProj, lcDiff, maskvalues = 0)
  sealsDiff <- classify(sealsDiff, seals7ToEsacci)

  lcProj <- terra::mask(lcBase, sealsDiff, maskvalues = NA, inverse = TRUE)
  lcProj <- cover(lcProj, sealsDiff, values = NA)

  dir.create(dirname(task$output), recursive = TRUE, showWarnings = FALSE)
  writeRaster(lcProj, task$output,
              overwrite = TRUE, filetype = "GTiff",
              datatype = "INT2U", gdal = c("COMPRESS=DEFLATE"))

  unlink(taskTmp, recursive = TRUE)
  message("Finished ", task$output)
  task$output
}

tasks <- buildTasks()
message("Found ", length(tasks), " SEALS land-cover files to process.")

nCores <- min(max(1, length(tasks)), 4)
runParallel <- identical(Sys.getenv("PREPARE_SEALS_PARALLEL", "0"), "1")
if (runParallel && nCores > 1) {
  cl <- parallel::makeCluster(nCores, type = "FORK",
                              outfile = file.path(dirOut, paste0("prepare_", as.integer(Sys.time()), ".out")))
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  foreach(task = tasks, .packages = "terra", .errorhandling = "stop") %dopar% {
    processTask(task)
  }
} else {
  message("Processing SEALS land-cover files serially.")
  lapply(tasks, processTask)
}

nloadFiles <- copyNLoadFiles()
readme <- writeReadme(tasks, nloadFiles)

unlink(dirTmp, recursive = TRUE)
message("Wrote README: ", readme)
message("Done.")
