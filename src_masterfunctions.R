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


# A wrapper for TGScover to go through all images and master.df and nine.pix tables
CreateCoverTable <- function(samplesize, googleimagefolder, cellnumbers, save.pixel.rgb = FALSE,
                             ninepixCSVname = "./nine_pix_df.csv", masterdfname = "./default_master", plot.9by9 = FALSE){
  
  ##########  First, make the table that values will be plugged into ###############################
  mx <- matrix(0, nrow = samplesize, ncol = 12)
  c.names <- c("Count", "pxl.num", "p.tree", "p.grass", "p.soil", "p.water", 
               "maj.truth", "SMA.tree", "SMA.grass", "SMA.soil", "SMA.water", "maj.SMA")
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
    image <- png.read.plot(image.name)
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
    master.df$nine.tree[t] <- mean(nine.cell.vector.tree)
    master.df$nine.grass[t] <- mean(nine.cell.vector.grass)
    master.df$nine.soil[t] <- mean(nine.cell.vector.soil)
  }
  
  # Save the 9-pixel df to csv
  write.csv(nine.pix.df, ninepixCSVname)
  write.csv(master.df, masterdfname)
  
  return(master.df)
}




addSMA <- function(SMAfolder, tiffname, treeband, grband, soilband, parkboundary,
                   groundtruthboundary, df, isHDR = FALSE){
  
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
  
  # Rescale values if HDR file
  if(isHDR){
    tr.max <- max(TREE[!is.na(TREE)])
    tr.min <- min(TREE[!is.na(TREE)])
    TREE[TREE > 0] <- TREE[TREE > 0]/tr.max*255
    TREE[TREE < 0] <- TREE[TREE < 0]/tr.min*255
    # TREE <- (TREE + abs(tr.min))/(tr.max - tr.min)*255
    tr.max <- max(GRASS[!is.na(GRASS)])
    tr.min <- min(GRASS[!is.na(GRASS)])
    GRASS[GRASS > 0] <- GRASS[GRASS > 0]/tr.max*255
    GRASS[GRASS < 0] <- GRASS[GRASS < 0]/tr.min*255
    # GRASS <- (GRASS + abs(tr.min))/(tr.max - tr.min)*255# Might want to do the same with negatives if this ends up overestimating tree cover
    tr.max <- max(SOIL[!is.na(SOIL)])
    tr.min <- min(SOIL[!is.na(SOIL)])
    SOIL[SOIL > 0] <- SOIL[SOIL > 0]/tr.max*255
    SOIL[SOIL < 0] <- SOIL[SOIL < 0]/tr.min*255
    # SOIL <- (SOIL + abs(tr.min))/(tr.max - tr.min)*255
  }
  
  plot(TREE, axes=F,box = F, main = "Trees SMA")
plot(GRASS, axes=F,box = F, main = "Grass SMA")
plot(SOIL, axes=F,box = F, main = "Soil SMA")
  
  TREE <- mask(TREE, groundtruthboundary) # if you try to do this in the clipTIF, it will result in a smaller raster and indexing won't work
  GRASS <- mask(GRASS, groundtruthboundary)
  SOIL <- mask(SOIL, groundtruthboundary)
  
  for(t in 1:samplesize){
    
    # Some cells are now NA because the classification excluded cover chances in some places, so switch to 0 
    if(is.na(TREE[cellnumbers[t]])){
      TREE[cellnumbers[t]] <- 0
    }
    
    if(is.na(GRASS[cellnumbers[t]])){
      GRASS[cellnumbers[t]] <- 0
    }
    
    if(is.na(SOIL[cellnumbers[t]])){
      SOIL[cellnumbers[t]] <- 0
    }
    
    # Get percent cover
    sma.tree <- TREE[cellnumbers[t]]/255  
    sma.grass <- GRASS[cellnumbers[t]]/255
    sma.soil <- SOIL[cellnumbers[t]]/255
    
    # Define major cover for SMA
    if(sma.tree > sma.grass & sma.tree > sma.soil){
      SMA.cov <- "Trees"
    } else if(sma.grass > sma.tree & sma.grass > sma.soil){
      SMA.cov <- "Grass"
    } else{SMA.cov <- "Soil"}
    
    # Plug in to df
    df$SMA.tree[t] <- sma.tree
    df$SMA.grass[t] <- sma.grass
    df$SMA.soil[t] <- sma.soil
    df$maj.SMA[t] <- SMA.cov 
    # Add the data to master.df
  }
  
  return(df)
  
}

addData <- function(variable, df, colname, samplesize, cellnumbers){
  for(t in 1:samplesize){
  
  
  
  
  
  
  
  
}

}




# image_prefix <- "LC81680602014034LGN00"
# image_folder <- "C:/Users/nagelki-4/Desktop/nagelki4/Grad School/Projects/EleTree Analysis/SMA/Landsat/Mpala/LC81680602014034-SC20160914165712/Original"

