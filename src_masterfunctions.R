#############################  LOAD LIBRARIES  ######################################################################
library(raster)
library(rgdal)
library(jpeg)
library(ggmap)
library(ggplot2)
library(RgoogleMaps)
library(png)
library(maps)
library(fields)
library(sp)



##################################  FUNCTIONS  #################################################################

# Read and display jpeg (requires jpeg package)
jpg.read.plot <- function(jpg.name){
  require(jpeg)
  image <- readJPEG(jpg.name)
  
  # get the dimensions of the image
  max.dim <- max(dim(image))
  x <- dim(image)[2]
  y <- dim(image)[1]
  
  plot(0:max.dim, type='n', axes = FALSE)
  rasterImage(image, 0,0,x,y)
  
  return(image)
}

# Read and displays png and moves it to an object
png.read.plot <- function(png.name){
  require(png)
  image <- readPNG(png.name)
  
  # get the dimensions of the image
  max.dim <- max(dim(image))
  x <- dim(image)[2]
  y <- dim(image)[1]
  
  #   plot(0:max.dim, type='n', axes = FALSE)
  #   rasterImage(image, 0,0,x,y)
  
  return(image)
}

# Calc MSAVI2 from a folder and Landsat file specication
MSAVI2_Landsat_calc <- function(path, prefix){
    
  # Identify mask and band file paths (we just want 4 and 5)
  landsat.bands <- list.files(path, pattern = prefix, full.names = TRUE)
  band4_path <- grep(pattern = "band4", landsat.bands, value = TRUE)
  band5_path <- grep(pattern = "band5", landsat.bands, value = TRUE)
  maskFile_path <- grep(pattern = "cfmask.tif", landsat.bands, value = TRUE) #two files have "cfmask" in them, so just use the order of files
  
  # Load the rasters
  mask <- raster(maskFile_path)
  band4 <- raster(band4_path)
  band5 <- raster(band5_path)
  
  # values below zero or >10000 are no good
  band4[band4 > 10000 | band4 < 0] <- NA
  band5[band5 > 10000 | band5 < 0] <- NA
  
  # Apply the mask to bands 4 and 5, setting cells with a mask value != 0 (clear) to NA (setting 0's to 1's and running this way saves about 4 sec, so 2:30ish overall)
  mask[mask != 0] <- NA
  mask <- mask + 1
  band5 <- (band5 * mask)
  band4 <- (band4 * mask)
  
  # Stack
  rstack <- stack(band5, band4)
  
  # Calc MSAVI2
  print("Starting MSAVI2")
  MSAVI2 <- calc(rstack, function(x) (2*x[,1] + 1 - sqrt((2*x[,1]+1)^2 - 8*(x[,1] - x[,2])))/2)
  #MSAVI2 <- (2*band5 + 1 - sqrt((2*band5+1)^2 - 8*(band5 - band4)))/2
  
  return(MSAVI2)
}

# Calc NDVI from a folder and Landsat file specication
NDVI_Landsat_calc <- function(path, prefix){
  
  # Identify mask and band file paths (we just want 4 and 5)
  landsat.bands <- list.files(path, pattern = prefix, full.names = TRUE)
  band4_path <- grep(pattern = "band4", landsat.bands, value = TRUE)
  band5_path <- grep(pattern = "band5", landsat.bands, value = TRUE)
  maskFile_path <- grep(pattern = "cfmask.tif", landsat.bands, value = TRUE) #two files have "cfmask" in them, so just use the order of files
  
  # Load the rasters
  mask <- raster(maskFile_path)
  band4 <- raster(band4_path)
  band5 <- raster(band5_path)
  
  # values below zero or >10000 are no good
  band4[band4 > 10000 | band4 < 0] <- NA
  band5[band5 > 10000 | band5 < 0] <- NA
  
  # Apply the mask to bands 4 and 5, setting cells with a mask value != 0 (clear) to NA (setting 0's to 1's and running this way saves about 4 sec, so 2:30ish overall)
  mask[mask != 0] <- NA
  mask <- mask + 1
  band5 <- (band5 * mask)
  band4 <- (band4 * mask)
  
  # Stack
  rstack <- stack(band5, band4)
  
  # Calc NDVI
  print("Starting NDVI")
  rstack <- stack(band5, band4)
  NDVI <- calc(rstack, function(x) (x[,1]-x[,2])/(x[,1]+x[,2])) # (this is twice as fast)
  #NDVI3 <- (band5 - band4)/(band5 + band4)
  
  return(NDVI)
}

# Stack all the bands in a Landsat image
stackLandsat <- function(path, prefix){
  # List the files
  landsat.bands <- list.files(path, pattern = prefix, full.names = TRUE)
  
  # Stack the bands
  for(i in 1:length(landsat.bands)){
    if(i == 1){
      band <- raster(landsat.bands[i])
      b.stack <- stack(band)
    }else{
      band <- raster(landsat.bands[i])
      b.stack <- stack(b.stack, band)
    }
  }
  return(b.stack)
}

# Plot a Landsat RGB (requires stackLandsat function)
rgbLandsat <- function(path = "", prefix = "", stackname = "", is.stacked = FALSE, r, b){
  if(nchar(path) > 0){
    stackname <- stackLandsat(path, prefix)
  }
  
  # Identify the red and blue bands
  rband <- grep(paste0("band", as.character(r)), names(stackname))
  bband <- grep(paste0("band", as.character(b)), names(stackname))
  # get just the RGB bands 
  rgbstack <- stack(stackname[[rband:bband]]) # assumes green is in the middle
  f.name <- names(rgbstack)[1]
  plotRGB(rgbstack, scale = 2^16, axes = TRUE, main = substr(f.name, 1, nchar(f.name)-9), stretch = "lin")
  
}

# Crop and mask a tif (can take a variable or a file name)
clipTIF <- function(path = "", tifname, bandnum = "", clipboundary){
  # Run these lines if a tif file has been given, otherwise skip
  if(nchar(path) > 0){
    ## Read in the tif
    tifname <- raster(paste0(path, tifname), band = bandnum)
  }
  
  # Crop them all down to the boundary extent
  ras <- crop(tifname, extent(clipboundary))
  
  # Clip to the boundary
  final.ras <- mask(ras, clipboundary)
  
  return(final.ras)
}

# Download a google image based on the boundary of a raster's pixel
getGoogleimage <- function(rastername, ras.pixel.num, destfilename, zoomlevel, typeofmap){
  
  # Extract that pixel
  pxl <- rasterFromCells(rastername, ras.pixel.num)
  # plot(pxl)
  
  ### Get the corner coordinates
  # List the 2 opposite corners
  c1 <- cbind(xmin(pxl), ymin(pxl))
  c2 <- cbind(xmax(pxl), ymax(pxl))
  cc <- rbind(c1, c2)
  # Find projection and make spatial points
  prj <- crs(pxl)
  utmcoor2<-SpatialPoints(cc, proj4string=prj)
  
  # Convert to lat long
  latlong <- spTransform(utmcoor2, CRS("+proj=longlat"))
  # Make extent objects for plugging into google map function
  ext2 <- extent(latlong)
  x.range <- c(ext2@xmin, ext2@xmax)
  y.range <- c(ext2@ymin, ext2@ymax)
  
  # DOWNLOAD THE IMAGE
  GetMap.bbox(x.range, y.range, destfile = destfilename, maptype = typeofmap, zoom = zoomlevel)
}

# image_prefix <- "LC81680602014034LGN00"
# image_folder <- "C:/Users/nagelki-4/Desktop/nagelki4/Grad School/Projects/EleTree Analysis/SMA/Landsat/Mpala/LC81680602014034-SC20160914165712/Original"

