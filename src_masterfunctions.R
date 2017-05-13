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
  
  plot(0:max.dim, type='n', axes = FALSE)
  rasterImage(image, 0,0,x,y)
  
  return(image)
}


# Calc MSAVI2 from a folder and Landsat file specication
MSAVI2_Landsat_calc <- function(path, prefix){ # should allow a shapefile in this to clip to first
  
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


# Outputs % tree, grass, soil and water (if included) cover when given rasters of red, green and blue bands
TGSWcover <- function(redband, greenband, blueband, save.class.pixel = FALSE, classify.water = FALSE,
                      imagesavefolder = "", pxlnumber){
  ## CLASSIFY
  # Sum the RGB values. Low values are dark areas which correspond to trees. 
  sum.im.pre <- greenband + redband + blueband
  
  # Calculate the Green-Red Vegetation Index
  grvd.pre <- (greenband - redband) / (greenband + redband)
  
  # Calculate summary stats for images to determine if it is mostly grass or soil
  #  This is done because images with a lot of grass or soil require different thresholds for classifying tree cover
  
  # Look at the brightness of each cell. Helps establish grass vs soil dominance
  bright <- redband + greenband + blueband 
  bright.area <- bright > 2
  bright.frac <- length(bright.area[bright.area == 1]) / length(bright.area)
  
  ### NOTE: in the bright and red.area calcs, in the original code, instead of using red.c etc, these were
  ######## defined by red, green and blue. If results are way off now, could change this back, but will need those
  ######## as inputs to the function as well.
  
  # Look at pixels with higher red values than green
  red.area <- redband > greenband
  red.frac <- length(red.area[red.area == 1]) / length(red.area) # 0.4391536
  
  # If including water, classify it first
  if(classify.water){
    # Classify water and include in classification if it takes up more than 5% of the image
    water <- sum.im.pre < 1.3 & grvd.pre < .03
    # plot(water, main = "Water")
    water.frac <- length(water[water == 1])/length(water[]) # calculate the fraction of water
    water.frac.thresh <- 0.05 # When water occupies what fraction of the scene will it be included in the classification? This was done because trees are sometimes classified as water, so only class water when it's fairly certain it's in the image. 
    if(water.frac > water.frac.thresh){
      water.mask <- water
      water.mask[water.mask == 1] <- 50
      holy.image.1 <- grvd.pre + water.mask
      holy.image.2 <- sum.im.pre + water.mask
      holy.image.1[holy.image.1 > 5] <- NA
      holy.image.2[holy.image.2 > 5] <- NA
      grvd.pre <- holy.image.1
      sum.im.pre <- holy.image.2
    }
  }
  
  # TREES    
  # This tests for areas with more red and brighter to help with classification bc trees under repped in soil dominated
  # and over repped in grass/tree dominated images
  if(red.frac >= .4 & bright.frac >= .5 ){
    trees <- sum.im.pre < 1.5
    tree.equation <- "sum.im.pre < 1.5"
  }else if(red.frac <= .01 & bright.frac <= .05){
    trees <- sum.im.pre < 1.1
    tree.equation <- "sum.im.pre < 1.1"
  }else{
    trees <- sum.im.pre < 1.3
    tree.equation <- "sum.im.pre < 1.3"
  }
  
  # Now remove the tree pixels from the image that will go into identifying soil and grass
  # create tree mask
  tree.mask <- trees
  tree.mask[tree.mask == 1] <- 50
  holy.image <- grvd.pre + tree.mask
  holy.image.3 <- sum.im.pre + tree.mask
  holy.image[holy.image > 5] <- NA
  holy.image.3[holy.image.3 > 5] <- NA
  grvd <- holy.image
  sum.im <- holy.image.3
  
  # GRASS
  # Split soil and grass using GRVD and brightness. Soil is nonveg and bright
  grass <-  sum.im < 2.1 & grvd > -0.05
  grass.equation <- "sum.im < 2.1 & grvd > -0.05"
  
  # SOIL 
  # Soil is the rest
  soil <- grass == 0
  soil.equation <- "remainder"
  
  # Combine into one
  # Give the different cover types different values
  grass.ad <- grass
  grass.ad[grass.ad == 1] <- 2
  grass.ad[is.na(grass.ad)] <- 0
  
  soil.ad <- soil # don't need to set soil to a different value. Just let it be 1
  soil.ad[is.na(soil.ad)] <- 0
  
  tree.ad <- trees
  tree.ad[tree.ad == 1] <- 3
  tree.ad[is.na(tree.ad)] <- 0
  
  par(mar=c(1,1,3,1))
  
  if(classify.water){
    if(water.frac > water.frac.thresh){
      water.ad <- water
      water.ad[water.ad == 1] <- 4
      water.ad[is.na(water.ad)] <- 0
      tot.img <- tree.ad + soil.ad + grass.ad + water.ad
      tot.img[tot.img == 0] <- NA
      
      # Save the class image, if wanted
      if(save.class.pixel){
        class.3030.name <- paste0("./", imagesavefolder, "/pixel_", pxlnumber, "_3030_class_new_redbrightthresh.png")
        # Save the image
        png(class.3030.name)
        plot(tot.img, axes=F,box = F, main = "4 = Water, 3 = Trees, 2 = Grass, 1 = Soil")
        dev.off()
      }
    } 
  }else{
    # Combine
    tot.img <- tree.ad + soil.ad + grass.ad
    tot.img[tot.img == 0] <- NA
    
    # Save the class image, if wanted
    if(save.class.pixel){
      class.3030.name <- paste0("./", imagesavefolder, "/pixel_", pxlnumber, "_3030_class_new_redbrightthresh.png")
      # Save the image
      png(class.3030.name)
      plot(tot.img, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")
      dev.off()
    }
  }
  
  # Calculate the % covers
  # Get the number of cells in each
  wtr.cells <- length(tot.img[tot.img == 4]) 
  tr.cells <- length(tot.img[tot.img == 3])
  gr.cells <- length(tot.img[tot.img == 2])
  sl.cells <- length(tot.img[tot.img == 1])
  tot.cells <- ncell(tot.img)
  
  # Compute percentages and what the majority cover is
  p.tree <- tr.cells / tot.cells
  p.grass <- gr.cells / tot.cells
  p.soil <- sl.cells / tot.cells
  p.water <- wtr.cells / tot.cells
  
  # Define majority cover for Google classification
  if(p.water > p.grass & p.water > p.soil & p.water > p.tree){
    maj.cov <- "Water"
  } else if(p.tree > p.grass & p.tree > p.soil & p.tree > p.water){
    maj.cov <- "Trees"
  } else if(p.grass > p.tree & p.grass > p.soil & p.grass > p.water){
    maj.cov <- "Grass"
  } else{maj.cov <- "Soil"}
  
  # Create the vector for output
  TGSWvec <- c(p.tree, p.grass, p.soil, p.water, maj.cov, tree.equation, grass.equation, soil.equation)
  names(TGSWvec) <- c("tree", "grass", "soil", "water", "majcov", "tree.equation", "grass.equation", "soil.equation")
  
  return(TGSWvec)
}

# Crops to 30 meters (google image). The function below (CreateCoverTable) does this itself, but could maybe incorporate
#  For use in download-create-record_GROUNDTRUTH.R
crop.google <- function(image.name, pixel.res = 30){
  require(png)
  # Now bring in the png that was just created
  image <- readPNG(auto.images[1])
  
  # Pull out the red green and blue bands
  red <- raster(image[,,1])
  green <- raster(image[,,2])
  blue <- raster(image[,,3])
  
  # Calculate how much to take off each end. The map is ~ 190x190 m
  image.indent <- (190 - pixel.res)/2
  # Crop the rasters
  one.side <- image.indent/190
  # Create the extent of the 30x30 pixel (the dimension can be changed in the variable section)
  cr.ext <- extent(c(one.side,1-one.side,one.side,1-one.side))
  
  # Crop
  red.c <- crop(red, cr.ext)
  green.c <- crop(green, cr.ext)
  blue.c <- crop(blue, cr.ext)
  
  # Stack and create the RGB image
  col.stack <- stack(red.c, green.c, blue.c)
  
  return(col.stack)
}

# A wrapper for TGScover to go through all images and master.df and nine.pix tables
CreateCoverTable <- function(samplesize, googleimagefolder, cellnumbers, Google.image.dim.final = 30, save.pixel.rgb = FALSE,
                             ninepixCSVname = "./nine_pix_df.csv", masterdfname = "./default_master", plot.9by9 = FALSE){
  
  require(png)
  
  ##########  First, make the table that values will be plugged into ###############################
  mx <- matrix(0, nrow = samplesize, ncol = 12)
  c.names <- c("Count", "pxl.num", "p.tree", "p.grass", "p.soil", "p.water", 
               "maj.truth", "SMA.tree", "SMA.grass", "SMA.soil", "SMA.water", 
               "maj.SMA")
  colnames(mx) <- c.names
  master.df <- as.data.frame(mx)
  
  # Make the 9-pixel table
  mx3 <- matrix(0, nrow = samplesize*9, ncol = 6)
  c3.names <- c("Count", "pxl.num", "number.of.9", "tree", "grass", "soil")
  colnames(mx3) <- c3.names
  nine.pix.df <- as.data.frame(mx3)
  
  
  ############# Define variables  ##################################################################
  # Start a 9 pixel count (will total the sample size * 9 in the end)
  count.nine.pix <- 1
  
  # Index list for the extents. Used to count through the 9 cells
  extent.list <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
  
  # Generate stats for indenting Google image to get the 9 different cells within it
  cell.size <- Google.image.dim.final/190 # 190m is the length of one side of the Google image
  half.pix <- cell.size / 2
  center <- .5 # Middle of the image. The reference point for indexing
  
  # Create the boundaries of the upper left cell. The 9-pixel loop will modify 
  # these as it counts through the cells (extent.list)
  bottom <- center + half.pix
  right <- center - half.pix
  top <- center + half.pix + cell.size 
  left <- center - half.pix - cell.size # might be able to simplify these. Write whole loop first
  
  
  ######## LOOP AND CLASSIFY  ##############################################################
  for(t in 1:samplesize){
    
    ###  Bring in pixel
    # Create the image name for that pixel
    image.name <- paste0("./", googleimagefolder, "/pixel_", cellnumbers[t], ".png")
    # Bring in the png
    image <- readPNG(image.name)
    # Pull out the red green and blue bands
    red <- raster(image[,,1])
    green <- raster(image[,,2])
    blue <- raster(image[,,3])
    
    # Create vectors to fill in the 9 values
    nine.cell.vector.tree <- c()
    nine.cell.vector.grass <- c()
    nine.cell.vector.soil <- c()
    
    
    #########  CLASSIFY 9-PIXEL and PLUG IN VALUES  and SAVE  ###########################################
    
    for(r in extent.list){ # This loop starts with the upper right of the nine cells, and reads down the first column, then second, then third
      
      # CREATE THE EXTENT
      # Calculate the extent of each cell. The map is ~ 190x190 m, ext(map) = 0, 1, 0, 1
      # First column
      if(r <= 2){
        bt <- bottom - cell.size*r
        rt <- right 
        tp <- top - cell.size*r
        lf <- left
      }else if(r <= 5){ # Second column
        bt <- bottom - cell.size*(r-3)
        rt <- right + cell.size 
        tp <- top - cell.size*(r-3)
        lf <- left + cell.size
      }else if(r <= 8){ # Third
        bt <- bottom - cell.size*(r-6)
        rt <- right + cell.size*2 
        tp <- top - cell.size*(r-6)
        lf <- left + cell.size*2
      }
      
      # Create extent object 
      cr.ext <- extent(c(lf, rt, bt, tp))
      
      if(plot.9by9){# Plot one of the layers to look at where the cells end up. Should be 3x3 grid
        plot(red)
        # Coerce extent to a SpatialPolygons object and add to plot
        cell.poly <- as(cr.ext, 'SpatialPolygons')
        plot(cell.poly, add = TRUE)
      }
      
      # Crop
      red.c <- crop(red, cr.ext)
      green.c <- crop(green, cr.ext)
      blue.c <- crop(blue, cr.ext)
      
      
      # Plot and save the RGB image of the center cell # not needed if they are already saved
      if(save.pixel.rgb & r == 4){
        col.stack <- stack(red.c, green.c, blue.c)
        rgb.3030.name <- paste0("./", googleimagefolder, "/pixel_", cellnumbers[t], "_3030_rgb.png")
        # Save the image
        png(rgb.3030.name)
        plotRGB(col.stack, scale = 1)
        dev.off()
      }
      
      
      
      # GET THE tree, grass and soil covers, along with the maj.cov and equations for the covers
      # Function written by Ryan. Function in src_masterfunctions.R (above)
      TGSWcovers <- TGSWcover(redband = red.c, blueband = blue.c, greenband = green.c, save.class.pixel = FALSE, classify.water = FALSE, 
                              imagesavefolder = googleimagefolder, pxlnumber = cellnumbers[t])
      
      
      
      # Record the covers
      nine.cell.vector.tree <- c(nine.cell.vector.tree, TGSWcovers[[1]])
      nine.cell.vector.grass <- c(nine.cell.vector.grass, TGSWcovers[[2]])
      nine.cell.vector.soil <- c(nine.cell.vector.soil, TGSWcovers[[3]])
      
      # Calc the cell number of the 9
      pix.num <- length(nine.cell.vector.tree)
      
      # Plug into dfs - master.df only when on the middle pixel (r == 4)
      if(r == 4){
        master.df$Count[t] <- t
        master.df$pxl.num[t] <- cellnumbers[t]
        master.df$p.tree[t] <- TGSWcovers[[1]]  #p.tree of the center of 9 pixels
        master.df$p.grass[t] <- TGSWcovers[[2]]
        master.df$p.soil[t] <- TGSWcovers[[3]]
        master.df$p.water[t] <- TGSWcovers[[4]]
        master.df$maj.truth[t] <- TGSWcovers[[5]]
        master.df$tree_equation[t] <- TGSWcovers[[6]]
        master.df$grass_equation[t] <- TGSWcovers[[7]]
        master.df$soil_equation[t] <- TGSWcovers[[8]]
      }
      
      # Fill the nine.pix
      nine.pix.df$Count[count.nine.pix] <- count.nine.pix
      nine.pix.df$pxl.num[count.nine.pix] <- cellnumbers[t]
      nine.pix.df$number.of.9[count.nine.pix] <- pix.num
      nine.pix.df$tree[count.nine.pix] <- TGSWcovers[[1]] # p.tree
      nine.pix.df$grass[count.nine.pix] <- TGSWcovers[[2]] # p.grass
      nine.pix.df$soil[count.nine.pix] <- TGSWcovers[[3]] # p.soil
      
      # Print
      print(count.nine.pix)
      
      # Add 1 to the count
      count.nine.pix <- count.nine.pix + 1
    }
    # Plug them in the averages for the 9 pixels into the df after all 9 are done
    master.df$nine.tree[t] <- mean(as.numeric(nine.cell.vector.tree))
    master.df$nine.grass[t] <- mean(as.numeric(nine.cell.vector.grass))
    master.df$nine.soil[t] <- mean(as.numeric(nine.cell.vector.soil))
  }
  
  # Save the 9-pixel df to csv
  write.csv(nine.pix.df, ninepixCSVname)
  write.csv(master.df, masterdfname)
  
  return(master.df)
}

# Add the SMA data to the master.df, and return that with the park level tree, grass, soil layers 
addSMA <- function(SMAfolder, tiffname, treeband, grband = "", soilband, parkboundary,
                   groundtruthboundary, df, isHDR = FALSE, is.unconstrained = FALSE){
  
  # get some stats
  samplesize <- nrow(df)
  cellnumbers <- df$pxl.num
  
  # Load and clip the SMA results to the park boundary
  # 2014
  TREE <- clipTIF(path = SMAfolder, tifname = tiffname, bandnum = treeband, clipboundary = parkboundary)
  GRASS <- clipTIF(path = SMAfolder, tifname = tiffname, bandnum = grband, clipboundary = parkboundary)
  SOIL <- clipTIF(path = SMAfolder, tifname = tiffname, bandnum = soilband, clipboundary = parkboundary)
  
  # 1987
  # TREE.87 <- clipTIF(path = SMA.folder, tifname = SMA.87, bandnum = tr.bandnum.87, clipboundary = mpala.boundary.simple)
  # GRASS.87 <- clipTIF(path = SMA.folder, tifname = SMA.87, bandnum = gr.bandnum.87, clipboundary = mpala.boundary.simple)
  # SOIL.87 <- clipTIF(path = SMA.folder, tifname = SMA.87, bandnum = so.bandnum.87, clipboundary = mpala.boundary.simple)
  
  # ENVI put zeros in for values that should be NA (They are values beyond the MESMA constraints)
  # A pure zero value is, I think, virtually impossible, so this shouldn't be overwriting any valid values
  TREE[TREE == 0] <- NA
  GRASS[GRASS == 0] <- NA
  SOIL[SOIL == 0] <- NA
  
  # Keep the raw values to plug into the df later
  raw.TREE <- TREE
  raw.GRASS <- GRASS
  raw.SOIL <- SOIL
  
  # Rescale values if HDR file
  if(isHDR){
    
    ## This is done to replicate how ENVI does it, but not saying it's the best way
    #### If the values are between 0-1, like they should be, this will do nothing to the HDR except scale to 255 (assuming max = 1 and min = 0)
    ### When they aren't 0-1, it creates scales the RGB255 between the min and max (that would be a problem if the min>0)
    ### See MASTER_Classifying Tree Cover.doc's HDR->Tiff section for more explanation 
    
    # Get the whole area that was classified, because those are the min and maxes needed for scaling to 255
    Laik.tree <- raster(paste0(SMAfolder, tiffname), band = treeband)
    Laik.grass <- raster(paste0(SMAfolder, tiffname), band = grband)
    Laik.soil <- raster(paste0(SMAfolder, tiffname), band = soilband)
    
    if(is.unconstrained == TRUE){
      # If it's unconstrained, scale it back to -1 to 2 range
      # Exclude all else
      Laik.tree[Laik.tree > 2] <- NA
      Laik.tree[Laik.tree < -1] <- NA
      Laik.grass[Laik.grass > 2] <- NA
      Laik.grass[Laik.grass < -1] <- NA
      Laik.soil[Laik.soil > 2] <- NA
      Laik.soil[Laik.soil < -1] <- NA
      TREE[TREE > 2] <- NA
      TREE[TREE < -1] <- NA
      GRASS[GRASS > 2] <- NA
      GRASS[GRASS < -1] <- NA
      SOIL[SOIL > 2] <- NA
      SOIL[SOIL < -1] <- NA
    }
    
    
    # Scale to 255 by shifting all to positive values and scaling to the max
    tr.min <- cellStats(Laik.tree, min)
    tr.max <- cellStats(Laik.tree, max)
    TREE.shift <- TREE - tr.min # brings the min to zero
    new.max <- tr.max - tr.min # find the new shifted max
    TREE <- (TREE.shift/new.max)*255
    # TREE <- (TREE - tr.min)/(tr.max - tr.min)*255 # This is the same as above, but simpler
    
    gr.min <- cellStats(Laik.grass, min)
    gr.max <- cellStats(Laik.grass, max)
    GRASS.shift <- GRASS - gr.min # brings the min to zero
    new.max <- gr.max - gr.min # find the new shifted max
    GRASS <- (GRASS.shift/new.max)*255
    
    
    so.min <- cellStats(Laik.soil, min)
    so.max <- cellStats(Laik.soil, max)
    SOIL.shift <- SOIL - so.min # brings the min to zero
    new.max <- so.max - so.min # find the new shifted max
    SOIL <- (SOIL.shift/new.max)*255
    
  }
  
  plot(TREE, axes=F,box = F, main = "Trees SMA")
  plot(GRASS, axes=F,box = F, main = "Grass SMA")
  plot(SOIL, axes=F,box = F, main = "Soil SMA")
  
  TREE.c <- mask(TREE, groundtruthboundary) # if you try to do this in the clipTIF, it will result in a smaller raster and indexing won't work
  GRASS.c <- mask(GRASS, groundtruthboundary)
  SOIL.c <- mask(SOIL, groundtruthboundary)
  rawT.c <- mask(raw.TREE, groundtruthboundary)
  rawG.c <- mask(raw.GRASS, groundtruthboundary)
  rawS.c <- mask(raw.SOIL, groundtruthboundary)
  
  for(t in 1:samplesize){
    
    # Some cells are now NA because the classification excluded cover chances
    # That's okay, but need a value for the major cover determination
    # So set % cover to zero for NA cells, else just calc the cover
    # if(is.na(TREE.c[cellnumbers[t]])){
    #   sma.tree <- 0
    # }else{sma.tree <- TREE.c[cellnumbers[t]]/255}
    # 
    # if(is.na(GRASS.c[cellnumbers[t]])){
    #   sma.grass <- 0
    # }else{sma.grass <- GRASS.c[cellnumbers[t]]/255}
    # 
    # if(is.na(SOIL.c[cellnumbers[t]])){
    #   sma.soil <- 0
    # }else{sma.soil <- SOIL.c[cellnumbers[t]]/255}
    # 
    # 
    
    
    if(TREE.c[cellnumbers[t]] >= 0.5){
      SMA.cov <- "Trees"
    }else{SMA.cov <- "Soil"}
    
    # # Define major cover for SMA
    # if(sma.tree > sma.grass & sma.tree > sma.soil){
    #   SMA.cov <- "Trees"
    # } else if(sma.grass > sma.tree & sma.grass > sma.soil){
    #   SMA.cov <- "Grass"
    # } else{SMA.cov <- "Soil"}
    # 
    # Plug in to df
    df$SMA.tree[t] <- TREE.c[cellnumbers[t]]#/255 # These stay in the original form because want NA cells to stay NA
    df$SMA.grass[t] <- GRASS.c[cellnumbers[t]]#/255 # That is to avoid them being plotted in the 1:1 plot
    df$SMA.soil[t] <- SOIL.c[cellnumbers[t]]#/255
    df$maj.SMA[t] <- SMA.cov 
    df$raw.SMA.tree[t] <- rawT.c[cellnumbers[t]]
    df$raw.SMA.grass[t] <- rawG.c[cellnumbers[t]]
    df$raw.SMA.soil[t] <- rawS.c[cellnumbers[t]]
  }
  
  r.list <- list(df, TREE, GRASS, SOIL)
  return(r.list)
  
}


addData <- function(tiffname, df, colname){
  # get some stats
  samplesize <- nrow(df)
  cellnumbers <- df$pxl.num
  
  # get the columb number from the df
  if(length(which(colnames(df) == colname)) < 1){ # if this is true, there is no column for the variable
    for(t in 1:samplesize){
      df$newcol[t] <- tiffname[cellnumbers[t]]
    }
    # Rename the new (last) column
    names(df)[ncol(df)] <- colname
  } else { # This assumes the column does exist then
    colnum <- which(colnames(df) == colname)
    
    for(t in 1:samplesize){
      df[t, colnum] <- tiffname[cellnumbers[t]]
    }
  }
  
  return(df)
  
}


confusionMatrix <- function(df, save.confusion.matrix = FALSE, save.matrix.as = "default_cnfusn_mtrx.csv", classify.water = FALSE){
  
  # get some stats
  sample.size <- nrow(df)
  
  if(classify.water){
    c.name <- c("Trees", "Grass", "Soil", "Water", "Total", "Percent Correct")
  }else{c.name <- c("Trees", "Grass", "Soil", "Total", "Percent Correct")}
  
  # make the table
  c.mtx <- matrix(0, nrow = length(c.name), ncol = length(c.name)) 
  c.mtx.df <- as.data.frame(c.mtx)
  colnames(c.mtx.df) <- c.name
  rownames(c.mtx.df) <- c.name
  
  # n <- 1
  for(n in 1:sample.size){
    truth <- df$maj.truth[n]
    pred <- df$maj.SMA[n]
    
    # Assign the coordinates for plugging into table
    if(truth == "Trees"){
      truth.cor <- 1
    } else if (truth == "Grass"){
      truth.cor <- 2
    } else if (truth == "Soil"){
      truth.cor <- 3
    }
    
    if(pred == "Trees"){
      pred.cor <- 1
    } else if (pred == "Grass"){
      pred.cor <- 2
    } else if (pred == "Soil"){
      pred.cor <- 3
    }
    
    # Now go the position in the table and add 1
    c.mtx.df[truth.cor, pred.cor] <- c.mtx.df[truth.cor, pred.cor] + 1
  }
  
  # Now sum the columns and calc the percent correct
  # First define how many columns/rows have numbers to add
  summing.values <- ncol(c.mtx.df) - 2
  tot.rowcol <- ncol(c.mtx.df) - 1
  per.rowcol <- ncol(c.mtx.df)
  tot.correct <- 0 # set up total correct variable to be adding values to from the diagonal 
  
  # Now calculate them
  for(i in 1:summing.values){
    c.mtx.df[i, tot.rowcol] <- sum(c.mtx.df[i, 1:summing.values]) # does the Total column
    c.mtx.df[i, per.rowcol] <- c.mtx.df[i, i] / c.mtx.df[i, tot.rowcol] # does the percent column
    c.mtx.df[tot.rowcol, i] <- sum(c.mtx.df[1:summing.values, i]) # does the Total row
    c.mtx.df[per.rowcol, i] <- c.mtx.df[i, i] / c.mtx.df[tot.rowcol, i] # does the percent column
    
    tot.correct <- tot.correct + c.mtx.df[i, i]
    # Fill in the total/total cell and the overall percent correct 
    if(i == summing.values){
      # Fill in the Total/Total combo cell
      c.mtx.df[tot.rowcol, tot.rowcol] <- sum(c.mtx.df[1:summing.values, tot.rowcol])
      # Last, fill in the total percent correct (lower right corner cell)
      c.mtx.df[per.rowcol, per.rowcol] <- tot.correct / c.mtx.df[tot.rowcol, tot.rowcol]
    }
  }
  
  if(save.confusion.matrix){
    write.csv(c.mtx.df, save.matrix.as)
  }
  
  return(c.mtx.df)
}


# Convert a variable name to text 
var2text <- function(variable){
  v.text <- deparse(substitute(variable))
  return(v.text)
}


# Plot 1 to 1 plot with a 1:1 line added and RMSE and Correlation in title
plot1to1 <- function(x, y, xlab = "", ylab = "", main = "", xlim = c(0,100), ylim = c(0,100), add_one2one_line = TRUE, 
                     save.plot = FALSE, save.plot.name = "rplot", col = "black", add.reg.line = "FALSE", hdr.table = "", hdr.count = "", 
                     hdr.namefor.table = ""){
  
  # Get rid of any rows that have an NA in either x or y
  x <- x[!is.na(y)]
  y <- y[!is.na(y)] # These have to be in this order so each removes the other's values based on the first one's NA values
  y <- y[!is.na(x)] # That is, you use the NAs in y to get rid of corresponding values/rows in x, then you get rid of the NA's in y
  x <- x[!is.na(x)] # Then the opposite is done with x's NAs getting rid of y values, then x getting rid of its own
  
  
  # Calculate Stats
  ##   RMSE
  rmse <- sqrt(mean((y - x)^2 , na.rm = TRUE))
  correl <- cor(x, y)
  reg <- lm(y ~ x) # this lm is just for reporting in the plot title and table
  
  # Plug values into table
  if(hdr.count > 0){
    
    hdr.table$hdr.name[hdr.count] <- hdr.namefor.table 
    hdr.table$RMSE[hdr.count] <- round(rmse, 2)
    hdr.table$adj.r.sq[hdr.count] <- round(summary(reg)[[9]], 2)
    hdr.table$slope[hdr.count] <- round(reg$coefficients[[2]], 2)
    hdr.table$'one-slope'[hdr.count] <- 1 - round(reg$coefficients[[2]], 2)
    hdr.table$'y-int'[hdr.count] <- round(reg$coefficients[[1]], 2)
    
  }
  
  
  
  # Do some controlling
  if(nchar(xlab) < 1) xlab <- var2text(x)
  if(nchar(ylab) < 1) ylab <- var2text(y)
  if(nchar(main) < 1) main <- paste(xlab, "vs.", ylab)
  
  # Plot the scaled 1 to 1
  plot(x*100, y*100, 
       main = paste0(main, " \nRMSE = ", round(rmse, 2)*100, "% "," R^2 = ", round(summary(reg)[[9]], 2) ,
                     " Slope = ", round(reg$coefficients[[2]], 2)), xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, col = col)
  if(add_one2one_line){
    abline(0,1) # Add the 1:1 line
  }
  if(add.reg.line){
    xnew <- x*100 # won't do regression on x*100 by itself
    ynew <- y*100
    modl <- lm(ynew ~ xnew) # calc regression
    abline(modl, lty = 2)
    # Add another line, shifted up to 1:1 line
    mdl2 <- modl
    mdl2$coefficients[[1]] <- 0
    abline(mdl2, lty = 2, col = "red")
  }
  
  if(save.plot){
    png(save.plot.name)
    plot(x*100, y*100, 
         main = paste0(main, " \nRMSE = ", round(rmse, 2)*100, "% "," R^2 = ", round(summary(reg)[[9]], 2) ,
                       " Slope = ", round(reg$coefficients[[2]], 2)), xlab = xlab, ylab = ylab, 
         xlim = xlim, ylim = ylim, col = col)
    if(add_one2one_line){
      abline(0,1) # Add the 1:1 line
    }
    if(add.reg.line){
      xnew <- x*100 # won't do regression on x*100 by itself
      ynew <- y*100
      modl <- lm(ynew ~ xnew) # calc regression
      abline(modl, lty = 2)
      # Add another line, shifted up to 1:1 line
      mdl2 <- modl
      mdl2$coefficients[[1]] <- 0
      abline(mdl2, lty = 2, col = "red")
    }
    dev.off()
  }
  if(hdr.count > 0){
    return(hdr.table)
  }
}


# Create's a transtion matrix from 2 stacks (one old and one new). Stack must be stack(Treecover, Grasscover, Soilcover) in that order
create.trans.mtx <- function(recent.TGS.cover.stack, yearofrecent = "", old.cover.TGS.cover.stack, yearofold = "", 
                             save.transition.matrix = FALSE, tran.mtx.name = "default.transition.matrix.csv"){
  
  # Change the names to be shorter
  r.stack <- recent.TGS.cover.stack
  o.stack <- old.cover.TGS.cover.stack
  
  # make the table
  t.mtx <- matrix(0, nrow = 3, ncol = 3)
  tran.mx <- as.data.frame(t.mtx)
  c.name <- c("Trees", "Grass", "Soil")
  colnames(tran.mx) <- c.name
  rownames(tran.mx) <- c.name
  
  
  
  # Plot RGB (Tree = Re)
  plotRGB(r.stack, scale = 255, main = paste0(yearofrecent, " \n", "Trees = Red, Grass = Green, Soil = Blue"))
  plotRGB(o.stack, scale = 255, main = paste0(yearofold, " \n", "Trees = Red, Grass = Green, Soil = Blue"))
  
  # Go to each cell in the two years and figure out dominant class
  # To start, set NA's to 0 so comparisons can be ran
  t87 <- o.stack[[1]]
  g87 <- o.stack[[2]]
  s87 <- o.stack[[3]]
  t14 <- r.stack[[1]]
  g14 <- r.stack[[2]]
  s14 <- r.stack[[3]]
  t87[is.na(t87)] <- 0
  g87[is.na(g87)] <- 0
  s87[is.na(s87)] <- 0
  t14[is.na(t14)] <- 0
  g14[is.na(g14)] <- 0
  s14[is.na(s14)] <- 0
  # n <- tran.sample[1]
  
  # Doing all the cells 
  for(n in 1:length(t87)){
    
    # Compute percentages and what the majority cover is
    tree_87 <- t87[n]
    grass_87 <- g87[n]
    soil_87 <- s87[n]
    tree_14 <- t14[n]
    grass_14 <- g14[n]
    soil_14 <- s14[n]
    
    # Define major cover types and skip to next loop if none is greater than the others (means they were all NA originally)
    # Define majority cover for 1987 SMA
    if(tree_87 > grass_87 & tree_87 > soil_87){
      cov.87 <- "Trees"
    } else if(grass_87 > tree_87 & grass_87 > soil_87){
      cov.87 <- "Grass"
    } else if(soil_87 > tree_87 & soil_87 > grass_87){
      cov.87 <- "Soil"
    } else{next}
    
    # Define major cover for 2014 SMA
    if(tree_14 > grass_14 & tree_14 > soil_14){
      cov.14 <- "Trees"
    } else if(grass_14 > tree_14 & grass_14 > soil_14){
      cov.14 <- "Grass"
    } else if(soil_14 > tree_14 & soil_14 > grass_14){
      cov.14 <- "Soil"
    } else{next}
    
    
    # Assign the table coordinates (row/column) for plugging into table
    if(cov.87 == "Trees"){
      coord.87 <- 1
    } else if (cov.87 == "Grass"){
      coord.87 <- 2
    } else if (cov.87 == "Soil"){
      coord.87 <- 3
    }
    
    if(cov.14 == "Trees"){
      coord.14 <- 1
    } else if (cov.14 == "Grass"){
      coord.14 <- 2
    } else if (cov.14 == "Soil"){
      coord.14 <- 3
    }
    
    # Now go to the position in the table and add 1
    tran.mx[coord.14, coord.87] <- tran.mx[coord.14, coord.87] + 1
  }
  
  # CREATE CHANGE MAPS
  # Subtract the three maps to get the change
  
  # Tree cover change
  tree.delta <- t14/255 - t87/255
  tree.delta[tree.delta[] == 0] <- NA
  plot(tree.delta, main = "Tree cover change")
  
  # Grass cover change
  grass.delta <- g14/255 - g87/255
  grass.delta[grass.delta[] == 0] <- NA
  plot(grass.delta, main = "Grass cover change")
  
  # Soil cover change
  soil.delta <- s14/255 - s87/255
  soil.delta[soil.delta[] == 0] <- NA
  plot(soil.delta, main = "Soil cover change")
  
  
  
  if(save.transition.matrix){
    write.csv(tran.mx, tran.mtx.name)
  }
  
  return(tran.mx)
  
  
  
  # 
  # # Create the entire list of cell numbers that have values
  # # Get the cell numbers
  # all.cell.val <- na.omit(mpala[])
  # 
  # all.cell.numbers <- c()
  # for(i in all.cell.val){
  #   sing.cell.new <- Which(mpala == i, cells= TRUE)
  #   all.cell.numbers <- c(all.cell.numbers, sing.cell.new)
  # }
  
  
  
  # # Check whether there are NA's within Mpala at all cover levels by summing them and if anything is still 0, then it was NA at all
  # t87 <- TREE.87
  # g87 <- GRASS.87
  # s87 <- SOIL.87
  # t87[is.na(t87)] <- 0
  # g87[is.na(g87)] <- 0
  # s87[is.na(s87)] <- 0
  # # Checked out. Would want to do this for the 14 year too, but thought not necessary
  # sum87 <- t87 + g87 + s87
  # plot(sum87)
  
  
}


# Output YearMonthDay
YMD <- function(){
  ymd <- Sys.time()
  ym <- substr(ymd, 1, 10)
  y <- gsub("-", "", ym)
  return(y)
}

# Move Landsat tar file and unzip in a new location. Keeps the original tar file in the new location
landsatUntar <- function(landsat.folder = "", landsat.file, dest.folder.name){
  # This is a simple function to unzip all the files in a folder (landsat.folder) in a new location
  # It can also take a single landsat file (landsat.file) and move it and unzip it
  
  # Load packages
  library(date)
  
  # if a folder name is provided, do the mass unzipping to the dest.folder
  if(nchar(landsat.folder) > 0){
    
    # Create list of files
    fileList <- list.files(landsat.folder, pattern = ".gz")
    numFiles <- length(fileList)
    
    
    for(i in 1:numFiles){
      # Use the date values to create Landsat_'year' folders and insert into those 
      # other folders with names in YYYYMMDD format within which to unzip the files
      # Print an update
      print(paste(i, "of", numFiles))
      
      # Create folder and file names
      # LandsatFolderName <- paste(j, "/Landsat_", year_num, sep = "")
      # PathRowFolderName <- paste(LandsatFolderName, "/Path", path_num, "_Row", row_num, sep = "")
      
      
      if (dir.exists(dest.folder.name) == FALSE) { # if folder doesn't exist, create it and move file
        dir.create(dest.folder.name)
        # Move file
        file.copy(fileList[i], dest.folder.name) # or do file.rename to cut and paste / moving before unzip might help speed
        untar(fileList[i], exdir = dest.folder.name)
      } else { # else just move the file to the appropriate folder
        file.copy(fileList[i], dest.folder.name)
        untar(fileList[i], exdir = dest.folder.name)
      }
      # Delete the tar.gz file after unzipping (saves space) (if untar didn't complete, then the file stays so can rerun code without missing files)
      file.remove(fileList[i])
    }
    
  }else{# else unzip the single file name to the dest.folder
    
    if (dir.exists(dest.folder.name) == FALSE) { # if folder doesn't exist, create it and move file
      dir.create(dest.folder.name)
      # Move file
      file.copy(landsat.file, dest.folder.name) # or do file.rename to cut and paste / moving before unzip might help speed
      untar(landsat.file, exdir = dest.folder.name)
    } else { # else just move the file to the appropriate folder
      file.copy(landsat.file, dest.folder.name)
      untar(landsat.file, exdir = dest.folder.name)
    }
    # Delete the tar.gz file after unzipping (saves space)
    file.remove(landsat.file)
  }
}


tick <- function(){
  return(Sys.time())
}

tock <- function(start.time){
  return(Sys.time() - start.time)
}

# Mosaics a list of Landsat bands, then stacks and clips to a boundary. 
# Can operate on a single file, too. It will then just stack it and clip it.
# The boundary needs to be a loaded shapefile 
# Will create a temp folder for the intermediate files
mosaic.stack.clip.mask.Landsat <- function(landsat.file.names, clip.boundary){
  
  # Required library
  library(raster)
  
  # Set the variable that deletes the intermediates to FALSE in case the next one up is just one image by itself
  delete.intermediates <- FALSE
  
  #########  Folder and file handling setup  ###################################################################
  # Create a temp folder in the working directory, if it doesn't already exist
  if(!file.exists("./temp")){
    dir.create("./temp") # all the intermeidate files will go here. For now, let's not delete them
  }
  
  
  #########  START  #######################################################################
  
  # What are the unique names?
  landsat.names <- substr(landsat.file.names, 1, 21)
  ls.unique <- unique(landsat.names)
  
  # Stitch if needed
  if(length(ls.unique) > 1){
    
    delete.intermediates <- FALSE # making this false for now. Later, when code is working well, could set it to TRUE
    
    # Merge (stitch) the seven bands and the cloudmask
    # Used merge instead of mosaic because mosaic blends the overlapping pixels. That would be and issue
    # in situations where there is cloud cover, especially for the cloud mask. Could have used mosaic
    # for the bands and merge for the cloud mask, but then if a band had clouds and was lower in the 
    # cloud mask merge, that mask wouldn't show up and the overlapping bands would have values between 
    # cloud and land surface values. So just merged them all so that the mask would always match what 
    # was showing up in the bands
    
    # For each band, rasterize and merge the images
    band.list <- c("band1", "band2", "band3", "band4", "band5", "band6", "band7", "pixel_qa")
    new.stitched.file.list <- ""
    # t <- 1
    for(t in 1:length(band.list)){
      # Get the names of the rasters with those bands
      same.path.n.bands <- grep(band.list[t], landsat.file.names, value = TRUE)
      
      # Create a text string for the different path and rows (will be part of the final name)
      pathrow.string <- ""
      
      ############  First, create the stitched file name and check to see if it already exists  #########################################
      for(n in 1:length(same.path.n.bands)){
        # update path row string
        pathrow.string <- paste0(pathrow.string, "_", substr(same.path.n.bands[n], 11, 16))
      }
      
      # This will be the name of the stitched file, so check if it already exists
      stitched.file.name <- paste0("./temp/", substr(same.path.n.bands[n], 1, 9), pathrow.string, substr(same.path.n.bands[n], 17, nchar(same.path.n.bands[n])))
      
      # If file exists, first add it to the name of tifs for the next section of code, then go to next in loop
      if(file.exists(stitched.file.name)){
        # Add tiff name to list
        new.stitched.file.list <- c(new.stitched.file.list, stitched.file.name)
        next
      }
      
      ########  For each raster, read it in and merge if it isn't the first raster  ###############################################
      for(n in 1:length(same.path.n.bands)){
        
        # update path row string
        pathrow.string <- paste0(pathrow.string, "_", substr(same.path.n.bands[n], 11, 16))
        
        if(n == 1){
          merge.raster <- raster(same.path.n.bands[n])
        }else{
          ras <- raster(same.path.n.bands[n])
          merge.raster <- merge(merge.raster, ras)
        }
      }
      
      # Write the stitched bands
      writeRaster(merge.raster, filename = stitched.file.name)
      
      # Add tiff name to list
      new.stitched.file.list <- c(new.stitched.file.list, stitched.file.name)
    }
    
    # Take off the first value in the list bc it will be ""
    landsat.file.names <- new.stitched.file.list[-1]
    
  }
  
  
  # For each Landsat image (there will only be one), stack it, clip it, make sure values are 0-10000 and apply the cloud mask. Then save as ENVI hdr
  
  # Get the 7 bands of the first image
  seven.band <- grep("band", landsat.file.names, value = TRUE)
  
  # Get the landsat prefix
  ls.prefix <- substr(seven.band[1], 1, nchar(seven.band[1]) - 13)
  
  # Get the cloud mask
  cloud.mask <- raster(paste0(ls.prefix, "_pixel_qa.tif"))
  
  # Set boundary to same crs as the landsat data
  site.boundary <- spTransform(clip.boundary, crs(cloud.mask)) 
  
  
  # Read in and stack the seven bands
  for(t in 1:length(seven.band)){
    ras <- raster(seven.band[t])
    if(t == 1){
      seven.stack <- ras
    }else{
      seven.stack <- stack(seven.stack, ras)
    }
  }
  
  # Crop to boundary
  clip.stack <- clipTIF(tifname = seven.stack, clipboundary = site.boundary)
  clip.mask <- clipTIF(tifname = cloud.mask, clipboundary = site.boundary)
  
  
  # For every band, set anything above 10000 or below 0 to NA
  for(t in 1:dim(clip.stack)[3]){
    # values below zero or >10000 are no good
    clip.stack[[t]][clip.stack[[t]] > 10000 | clip.stack[[t]] < 0] <- NA
    
  }
  
  
  # new qa band for clouds. description here: https://landsat.usgs.gov/collectionqualityband
  # band is called pixel_qa and 66 = clear (no water or clouds or snow)
  
  # Apply the cloud mask
  clip.mask[clip.mask != 66] <- NA
  mask <- clip.mask - 65 # values of 66 mean clear skies, so bring that to 1 so that when multiplied, those places stay the same
  # Apply mask to each layer of stack
  for(n in 1:dim(clip.stack)[3]){
    clip.stack[[n]] <- (clip.stack[[n]] * mask)   
  }
  
  
  # If images were stitched, this will delete the intermediates
  if(delete.intermediates == TRUE){
    file.remove(landsat.file.names)
  }
  
  return(clip.stack)
  
}


# This on puts them in a brick, instead of a stack, then clips, masks, sets values, and merges last
# No temp folder created
mosaic.stack.clip.mask.Landsat2 <- function(landsat.file.names, clip.boundary){
  
  # Required library
  library(raster)
  
  # Set the variable that deletes the intermediates to FALSE in case the next one up is just one image by itself
  delete.intermediates <- FALSE
  
  #########  Folder and file handling setup  ###################################################################
  # Create a temp folder in the working directory, if it doesn't already exist
  if(!file.exists("./temp")){
    dir.create("./temp") # all the intermeidate files will go here. For now, let's not delete them
  }
  
  
  #########  START  #######################################################################
  
  # What are the unique names?
  landsat.names <- substr(landsat.file.names, 1, 21)
  ls.unique <- unique(landsat.names)
  
  
  ### Go through each landsat image and merge at the end
  x <- 1
  for(x in 1:length(ls.unique)){
    
    ###############  BRICK  ##########################################################################################
    # For each Landsat image (there will only be one), brick it, clip it, make sure values are 0-10000 and apply the cloud mask. Then save as ENVI hdr
    
    # Get the 7 bands of the first image
    s.band <- grep(ls.unique[x], landsat.file.names, value = TRUE)
    seven.band <- grep("band", s.band, value = TRUE)
    
    
    # Read in and stack the seven bands
    for(t in 1:length(seven.band)){
      ras <- raster(seven.band[t])
      if(t == 1){
        s.stack <- ras
      }else{
        s.stack <- stack(s.stack, ras)
      }
    }
    
    # Conver to a brick
    s.stack <-  brick(s.stack)
    
    ###########  Cloud Mask and Reproject Boundary  ############################################################################# 
    
    # Get the cloud mask
    cloud.mask <- raster(grep("_pixel_qa.tif", s.band, value = TRUE))
    
    # Set boundary to same crs as the landsat data
    site.boundary <- spTransform(clip.boundary, crs(cloud.mask)) 
    
    
    ############  CLIP  ##########################################################################################
    # Crop to boundary
    s.stack <- clipTIF(tifname = s.stack, clipboundary = site.boundary)
    clip.mask <- clipTIF(tifname = cloud.mask, clipboundary = site.boundary)
    
    
    ###########  MASK  ################################################################################################
    clip.mask[clip.mask != 66] <- NA
    mask <- clip.mask - 65 # values of 66 mean clear skies, so bring that to 1 so that when multiplied, those places stay the same
    # Apply mask to the brick
    s.stack <- s.stack * mask
    
    
    # rs2 <- calc(clip.stack, fun=function(x){x * mask})
    # 
    
    #########  CHANGE HIGH and LOW VALUES to NA  ################################################################################
    # For every band, set anything above 10000 or below 0 to NA
    for(t in 1:dim(s.stack)[3]){
      # values below zero or >10000 are no good
      s.stack[[t]][s.stack[[t]] > 10000 | s.stack[[t]] < 0] <- NA
    }
    
    
    
    
    ##########  MERGE  ####################################################################################################
    # Stitch if needed, if not, end function
    
    if(length(ls.unique) == 1){
      print("Single image (no merging).")
      return(s.stack) 
      break
    }
    
    
    if(length(ls.unique) > 1){
      
      if(x == 1){
        first.brick <- s.stack
        brick.list <- list(first.brick)
      }else if(x == 2){
        second.brick <- s.stack
        brick.list <- list(first.brick, second.brick)
      }else if(x == 3){
        third.brick <- s.stack
        brick.list <- list(first.brick, second.brick, third.brick)
      }else if(x == 4){
        fourth.brick <- s.stack
        brick.list <- list(first.brick, second.brick, third.brick, fourth.brick)
      }else if(x == 5){
        fifth.brick <- s.stack
        brick.list <- list(first.brick, second.brick, third.brick, fourth.brick, fifth.brick)
      }
    }
    # if the merge went here, then it would merge after each iteration of ls.unique. So if there were 3 images, it would merge twice instead of once
  }
  
  # Just in case it doesn't break out of the function above with ls.unique == 1, this is written as an if statement
  if(length(ls.unique) > 1){
    merge.start <- tick()
    s.stack <- do.call(merge, brick.list) 
    print("Merge time:")
    print(tock(merge.start))
    
    return(s.stack)
  }
}


# From a list of years, this finds what year is in the name of something
findYear <- function(possible.years, files.to.search){
  # Start a running list
  temp <- c()
  for(i in files.to.search){
    for(x in possible.years){
      if(grepl(x, i)){
        temp <- c(temp, x)
      }
    }
  }
  return(temp)
}      

