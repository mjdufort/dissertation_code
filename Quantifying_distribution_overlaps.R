## This includes my scripts to quantify geographic range sizes and range overlaps of Picidae species in R, using the rgeos package.
## it also includes scripts to generate maps of species richness

library(sp)
library(rgeos)
library(rgdal)
library(shapefiles)
library(maptools)
library(maps)
library(rangeMapper)
library(latticeExtra)


### generate some useful functions (these may need tweaking for the migratory species)

## the function getAreas() determines the area of a list of spatial polygons
# as input, it takes spgeom.list (a list of sp objects)
# it returns a vector with the total area of each element of the list
getAreas <- function(spgeom.list, quiet=TRUE) {
  require(rgeos)
  areas <- rep(NA, length(spgeom.list))  # intialize a vector of areas
  for (i in 1:length(spgeom.list)) {
    if (!quiet) cat("Calculating area of item ", names(spgeom.list)[i], ".\n", sep="")
    areas[i] <- gArea(spgeom.list[[i]])  # calculate area of the current geometry
  }
  names(areas) <- names(spgeom.list)  # assign names of the elements in the list to the elements of the vector of areas
  return(areas)
}

## the function calculateOverlaps calculates the pairwise overlap of all range maps
# as input, it takes geom.ranges (a list of spatial polygons in sp format)
# it returns a matrix or list of matrices of pairwise overlaps (species in row x overlaps with species in column y)
# the optional argument mode allows to select absolute overlap (mode="absolute"), overlap scaled by the range of the species in the given row (mode="scaled"), or both (mode="both")
calculateOverlaps <- function(geom.ranges, mode="both", quiet=TRUE, alpha_order=TRUE) {
  require(rgeos)
  if (!(mode %in% c("scaled", "absolute", "both"))) {
    print(paste("Mode ",mode," not recognized", sep=""))
    return()
  }
  nranges <- length(geom.ranges)  # calculate the number of geometries
  taxnames <- names(geom.ranges)  # extract the names of the geometries
  overlaps <- matrix(nrow=nranges, ncol=nranges, dimnames=list(taxnames, taxnames))  # initialize an empty matrix to store the pairwise overlaps
  for (i in 1:nranges) {  # loop over geometries
    if (!quiet) cat("\nCalculating overlaps for ", taxnames[i], ".\n\n", sep="")
    area.tmp <- gArea(geom.ranges[[i]])  # calculate the area of the current geometry
    for (j in 1:i) {  # loop over geometries (again)
      if (!quiet) cat("Calculating overlap with ", taxnames[j], ".\n", sep="")
      if (i==j) {
        overlaps[i,j] <- area.tmp  # calculate overlap with self (geometry area)
      } else if (!gOverlaps(geom.ranges[[i]], geom.ranges[[j]])) {
        overlaps[i,j] <- 0  # if geometries don't overlap, set overlap to 0 (gIntersection fails in this case)
      } else {
        geom.tmp <- gIntersection(geom.ranges[[i]], geom.ranges[[j]])  # calculate the intersection of the two geometries
        if (!gIsValid(geom.tmp)) { 
          overlaps[i,j] <- NA  # if the intersection is not valid, set the overlap to NA
        } else overlaps[i,j] <- gArea(geom.tmp)  # if the intersection is valid, set the overlap to the area of the intersection
      }
      overlaps[j,i] <- overlaps[i,j]  # fill in the other half of the matrix
    }
  }
  
  if (alpha_order) overlaps <- overlaps[order(rownames(overlaps)), order(colnames(overlaps))]  # put the geometries in alphabetical order (if specified)
  
  # return the appropriate values
  if (mode %in% c("scaled","both")) overlaps.scaled <- overlaps / diag(overlaps)
  switch(mode,
         absolute = return(overlaps),
         scaled = return(overlaps.scaled),
         both = return(list(overlaps=overlaps, overlaps.scaled=overlaps.scaled))
         )
}

## the function calculateOverlapsMigratory() calculates overlaps between migratory and nonmigratory species
# as input, it takes geom.ranges.nonmigratory (a list of geometries for nonmigratory taxa), geom.ranges.migratory (a list of summer and winter geometries for migratory taxa), mode (option to calculate absolute overlaps, scaled overlaps, or both), and summer_fraction (the proportion of the overlap from the summer range; winter proportion is 1 - summer_fraction)
# it assumes that all taxa have consistent use of summer (e.g. austral winter would be considered summer, so that the timing is the same); this is important because migratory species can only overlap each other in their shared summer ranges (summer only compared to summer) or shared winter ranges
calculateOverlapsMigratory <- function(geom.ranges.nonmigratory, geom.ranges.migratory, mode="both", summer_fraction=0.5, quiet=TRUE, alpha_order=TRUE) {
  library(rgeos)
  
  if (!(mode %in% c("scaled", "absolute", "both"))) {  # check for invalid modes
    cat("Mode ", mode ," not recognized.\n", sep="")
    return()
  }

  # get number and names of geometries for nonmigratory species
  nranges.nonmigratory <- length(geom.ranges.nonmigratory)
  taxnames.nonmigratory <- names(geom.ranges.nonmigratory)
  
  # get number and names of geometries for migratory species
  nranges.migratory <- length(geom.ranges.migratory) / 2
  taxnames.migratory <- unique(gsub("_summer", "", gsub("_winter", "", names(geom.ranges.migratory))))  # pull out the taxon names for the migratory taxa
  taxnames.migratory
  
  # check list of ranges for incorrect numbers or taxa in both migratory and nonmigratory
  for (taxname in taxnames.migratory) {
    if (length(grep(taxname, names(geom.ranges.migratory))) != 2) {  # check if the number of ranges is consistent with assumptions (each migratory taxon should have 2!)
      cat("Migratory ranges are not valid. Function librarys summer and winter range for each migratory taxon.\n")
      return()
    }
    if (length(grep(taxname, names(geom.ranges.nonmigratory))) > 0) {  # check if the migratory taxon is also listed as a nonmigratory taxon
      cat("Error: ", taxname, " found in both nonmigratory and migratory ranges.\n", sep="")
      return()
    }
  }
  
  overlaps <- matrix(nrow=(nranges.nonmigratory + nranges.migratory), ncol=(nranges.nonmigratory + nranges.migratory), dimnames=list(c(taxnames.nonmigratory, taxnames.migratory), c(taxnames.nonmigratory, taxnames.migratory)))  # initialize empty matrix of overlaps
  overlaps.scaled <- overlaps  # create a copy to store scaled overlaps (because they need to be calculated as we go for migratory species)
  
  # output number of nonmigratory and migratory ranges
  cat("Nonmigratory ranges: ", nranges.nonmigratory, "\n", sep="")
  cat("Migratory ranges: ", nranges.migratory, "\n", sep="")
  
  # calculate areas for all geometries (nonmigratory taxa and summer and winter ranges for migratory taxa)
  areas <- as.numeric(rep(NA, length(geom.ranges.nonmigratory) + length(geom.ranges.migratory)))
  for (i in 1:nranges.nonmigratory) {
    areas[i] <- gArea(geom.ranges.nonmigratory[[i]])
  }
  for (j in 1:length(geom.ranges.migratory)) {
    areas[j+nranges.nonmigratory] <- gArea(geom.ranges.migratory[[j]])
  }
  names(areas) <- c(names(geom.ranges.nonmigratory), names(geom.ranges.migratory))
  
  # loop over the nonmigratory taxon x nonmigratory taxon pairs
  cat("\nCalculating overlaps among nonmigratory taxa.\n")
  for (i in 1:nranges.nonmigratory) {  # loop over first taxon
    if (!quiet) cat("\nCalculating overlaps for ", taxnames.nonmigratory[i], " with nonmigratory taxa.\n\n", sep="")
    area.tmp <- areas[[taxnames.nonmigratory[i]]]  # retrieve area for current taxon 2
    for (j in 1:i) {  # loop over second taxon (up through and including first taxon)
      area.tmp.target <- areas[[taxnames.nonmigratory[j]]]  # retrieve area for current taxon 2
      if (!quiet) cat("Calculating overlap with ", taxnames.nonmigratory[j], ".\n", sep="")
      if (i==j) {
        overlaps[i,j] <- area.tmp  # if calculating overlap with self, use range area
      } else if (!gOverlaps(geom.ranges.nonmigratory[[i]], geom.ranges.nonmigratory[[j]])) {
        overlaps[i,j] <- 0  # if ranges don't overlap, set overap to 0 (gIntersection fails in this case)
      } else {
        geom.tmp <- gIntersection(geom.ranges.nonmigratory[[i]], geom.ranges.nonmigratory[[j]])  # calculate intersection of the two ranges
        if (!gIsValid(geom.tmp)) {  # check if intersection is a valid geometry
          overlaps[i,j] <- NA  # if intersection is not a valid geometry, set overlap to 0
        } else overlaps[i,j] <- gArea(geom.tmp)  # if intersection is valid, set overlap to area of the intersection
      }
      
      overlaps[j,i] <- overlaps[i,j]  # fill in the other half of the matrix
      
      overlaps.scaled[i,j] <- overlaps[i,j] / area.tmp  # calculate scaled overlap
      overlaps.scaled[j,i] <- overlaps[i,j] / area.tmp.target  # fill in the other half of the matrix for scaled overlap (scaled by taxon 2)
    }
  }
  
  # loop over the migratory taxon x nonmigratory taxon pairs
  cat("\nCalculating overlaps between migratory and nonmigratory taxa.\n")
  for (i in 1:nranges.migratory) {  # loop over first taxon
    if (!quiet) cat("\nCalculating overlaps for ", taxnames.migratory[i], " with nonmigratory taxa.\n\n", sep="")
    area.tmp.summer <- areas[[paste(taxnames.migratory[i], "_summer", sep="")]]  # retrieve size of summer range
    area.tmp.winter <- areas[[paste(taxnames.migratory[i], "_winter", sep="")]]  # retrieve size of winter range
    
    for (j in 1:nranges.nonmigratory) {  # loop over second taxon
      if (!quiet) cat("Calculating overlap with ", taxnames.nonmigratory[j], ".\n", sep="")
      area.tmp.target <- areas[[taxnames.nonmigratory[j]]]  # retrieve size of taxon 2 range
      
      # calculate summer range overlap
      if (!gOverlaps(geom.ranges.migratory[[paste(taxnames.migratory[i], "_summer", sep="")]], geom.ranges.nonmigratory[[j]])) {  # check if summer range overlaps range of nonmigratory species
        overlap.tmp.summer <- 0  # if summer range doesn't overlap, set summer overlap to 0
      } else {
        geom.tmp.summer <- gIntersection(geom.ranges.migratory[[paste(taxnames.migratory[i], "_summer", sep="")]], geom.ranges.nonmigratory[[j]])  # calculate intersection of summer range and taxon 2 range
        if (!gIsValid(geom.tmp.summer)) {  # check if summer overlap is a valid geometry
          overlap.tmp.summer <- NA  # if intersection is not valid, set overlap to NA
        } else {
          overlap.tmp.summer <- gArea(geom.tmp.summer)  # if intersection is valid, set overlap to area of intersection
        }
      }
      
      # calculate winter range overlap
      if (!gOverlaps(geom.ranges.migratory[[paste(taxnames.migratory[i], "_winter", sep="")]], geom.ranges.nonmigratory[[j]])) {  # check if winter range overlaps range of nonmigratory species
        overlap.tmp.winter <- 0  # if winter range doesn't overlap, set winter overlap to 0
      } else {
        geom.tmp.winter <- gIntersection(geom.ranges.migratory[[paste(taxnames.migratory[i], "_winter", sep="")]], geom.ranges.nonmigratory[[j]])  # calculate intersection of winter range and taxon 2 range
        if (!gIsValid(geom.tmp.winter)) {  # check if winter overlap is a valid geometry
          overlap.tmp.winter <- NA  # if intersection is not valid, set overlap to NA
        } else {
          overlap.tmp.winter <- gArea(geom.tmp.winter)  # if intersection is valid, set overlap to area of intersection
        }
      }

      overlaps[i+nranges.nonmigratory, j] <- (overlap.tmp.summer * summer_fraction) + (overlap.tmp.winter * (1 - summer_fraction))  # store overlap
      overlaps[j, i+nranges.nonmigratory] <- overlaps[i+nranges.nonmigratory, j]  # fill in the other half of the matrix
      
      overlaps.scaled[i+nranges.nonmigratory, j] <- ((overlap.tmp.summer / area.tmp.summer) * summer_fraction) +   ((overlap.tmp.winter / area.tmp.winter) * (1 - summer_fraction))  # calculate overlap scaled by taxon 1; overlap in each season is scaled by the range size for that season
      overlaps.scaled[j, i+nranges.nonmigratory] <- overlaps[i+nranges.nonmigratory, j] / area.tmp.target  # fill in the other half of the matrix for scaled overlap (scaled by taxon 2 range)
      
    }
  }
  
  # loop over the migratory taxon x migratory taxon pairs
  cat("\nCalculating overlaps among migratory taxa.\n")
  for (i in 1:nranges.migratory) {  # loop over first taxon
    if (!quiet) cat("\nCalculating overlaps for ", taxnames.migratory[i], " with migratory taxa.\n\n", sep="")
    area.tmp.summer <- areas[[paste(taxnames.migratory[i], "_summer", sep="")]]  # retrieve size of summer range
    area.tmp.winter <- areas[[paste(taxnames.migratory[i], "_winter", sep="")]]  # retrieve size of winter range
    
    for (j in 1:i) {  # loop over second taxon
      if (!quiet) cat("Calculating overlap with ", taxnames.migratory[j], ".\n", sep="")
      area.tmp.target.summer <- areas[[paste(taxnames.migratory[j], "_summer", sep="")]]  # retrieve size of summer range of taxon 2
      area.tmp.target.winter <- areas[[paste(taxnames.migratory[j], "_winter", sep="")]]  # retrieve size of winter range of taxon 2
      
      if (i==j) {  # check if this case is overlap with self
        overlaps[i+nranges.nonmigratory, j+nranges.nonmigratory] <- (area.tmp.summer * summer_fraction) + (area.tmp.winter * (1-summer_fraction))  # if calculating overlap with self, set overlap to weighted average of summer and winter range size
        overlaps.scaled[i+nranges.nonmigratory, j+nranges.nonmigratory] <- 1  # if calculating overlap with self, set scaled overlap to 1
      } else {
        
        # calculate summer range overlaps
        if (!gOverlaps(geom.ranges.migratory[[paste(taxnames.migratory[i], "_summer", sep="")]], geom.ranges.migratory[[paste(taxnames.migratory[j], "_summer", sep="")]])) {  # check if summer ranges overlap
          overlap.tmp.summer <- 0
        } else {
          geom.tmp.summer <- gIntersection(geom.ranges.migratory[[paste(taxnames.migratory[i], "_summer", sep="")]], geom.ranges.migratory[[paste(taxnames.migratory[j], "_summer", sep="")]])
          if (!gIsValid(geom.tmp.summer)) {  # check if summer overlap is a valid geometry
            overlap.tmp.summer <- NA
          } else {
            overlap.tmp.summer <- gArea(geom.tmp.summer)
          }
        }
        
        # calculate winter range overlaps
        if (!gOverlaps(geom.ranges.migratory[[paste(taxnames.migratory[i], "_winter", sep="")]], geom.ranges.migratory[[paste(taxnames.migratory[j], "_winter", sep="")]])) {  # check if winter ranges overlap
          overlap.tmp.winter <- 0
        } else {
          geom.tmp.winter <- gIntersection(geom.ranges.migratory[[paste(taxnames.migratory[i], "_winter", sep="")]], geom.ranges.migratory[[paste(taxnames.migratory[j], "_winter", sep="")]])
          if (!gIsValid(geom.tmp.winter)) {  # check if winter overlap is a valid geometry
            overlap.tmp.winter <- NA
          } else {
            overlap.tmp.winter <- gArea(geom.tmp.winter)
          }
        }
        
        overlaps[i+nranges.nonmigratory, j+nranges.nonmigratory] <- (overlap.tmp.summer * summer_fraction) + (overlap.tmp.winter * (1 - summer_fraction))  # calculate overlaps as weighted average of summer and winter overlaps
        overlaps[j+nranges.nonmigratory, i+nranges.nonmigratory] <- overlaps[i+nranges.nonmigratory, j+nranges.nonmigratory]  # fill in the other half of the matrix
        overlaps.scaled[i+nranges.nonmigratory, j+nranges.nonmigratory] <- ((overlap.tmp.summer / area.tmp.summer) * summer_fraction) +   ((overlap.tmp.winter / area.tmp.winter) * (1 - summer_fraction))  # calculate scaled overlap, scaled by taxon 1 range
        overlaps.scaled[j+nranges.nonmigratory, i+nranges.nonmigratory] <- ((overlap.tmp.summer / area.tmp.target.summer) * summer_fraction) +   ((overlap.tmp.winter / area.tmp.target.winter) * (1 - summer_fraction))  # calculate scaled overlap for the other half of the matrix, scaled by taxon 2 range
      }
      
    }
  }
  
  # put matrix elements in alphabetical order (if specified)
  if (alpha_order) {
    overlaps <- overlaps[order(rownames(overlaps)), order(colnames(overlaps))]
    overlaps.scaled <- overlaps.scaled[order(rownames(overlaps.scaled)), order(colnames(overlaps.scaled))]
  }
  
  # return the appropriate values
  switch(mode,
         absolute = return(overlaps),
         scaled = return(overlaps.scaled),
         both = return(list(overlaps=overlaps, overlaps.scaled=overlaps.scaled))
  )
}


## the function checkValidity() checks the validity of a list of spatial polygons
# as input, it takes spgeom.list (a list of geometries as sp objects)
# it returns a vector of TRUE/FALSE
checkValidity <- function(spgeom.list, quiet=TRUE) {
  library(rgeos)
  validity <- rep(NA, length(spgeom.list))  # initialize a vector to store validity
  for (i in 1:length(spgeom.list)) {  # loop over geometries in list
    if (!quiet) cat("Checking validity for ", names(spgeom.list)[i], ".\n", sep="")
    validity[i] <- gIsValid(spgeom.list[[i]])  # test if current geometry is valid
  }
  names(validity) <- names(spgeom.list)  # assign geometry names to vector elements
  return(validity)
}


### read in the shapefiles and do some basic processing
## the process:
# 1: read in a list of the filenames
# 2: read in the shapefiles
# 3: set the CRS for the shapefiles
# 4: transform the shapefiles to equal area
# 5: simplify the geometries using gUnaryUnion; this DRAMATICALLY speeds things up, because it eliminates the need for the program to loop over many separate features (hundreds in some cases); also, need to do gUnaryUnion BEFORE gBuffer, or gBuffer fills in holes in the geometries (which we don't want!!!)
# 6: if any of them fail because of geometry problems, add a 0-width buffer using gUnaryUnion(gBuffer(problem_geometry, width=0)), or try to fix in an external GIS program
# 7: make any invalid geometries valid with gBuffer(, width=0); do this after gUnaryUnion because that helps prevent large changes in the geometry


## read in the geometries

taxon_list.BirdLife <- as.character(read.table(file="shp_list_BirdLife_all_year.txt")[,1])  # read in a list of the filenames (created externally)

shp_list.BirdLife <- list()  # initialize empty list
shp_list.BirdLife[["all_year"]] <- as.character(read.table(file="shp_list_BirdLife_all_year.txt")[,1])  # read list of filenames treating all species as nonmigratory
shp_list.BirdLife[["migratory"]] <- as.character(read.table(file="shp_list_BirdLife_migratory.txt")[,1])  # read list of filenames treating all species as nonmigratory

ranges.shp.BirdLife <- list()  # initialize a list to store the polygons in
## retrieve the polygons from shape files
for (i in names(shp_list.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  ranges.shp.BirdLife[[i]] <- list()
  for (j in 1:length(shp_list.BirdLife[[i]])) {  # loop over taxon
    ranges.shp.BirdLife[[i]][j] <- readShapePoly(shp_list.BirdLife[[i]][j], delete_null_obj=TRUE)
  }
  names(ranges.shp.BirdLife[[i]]) <- shp_list.BirdLife[[i]]
}
rm(i,j)

## set CRS for all of the polygons
for (i in names(ranges.shp.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  for (j in names(ranges.shp.BirdLife[[i]])) {  # loop over taxa
    proj4string(ranges.shp.BirdLife[[i]][[j]]) <- CRS("+init=epsg:4326")
  }
}
rm(i,j)

## check validity of the polygons
validity.ranges.shp.BirdLife <- list()
for (i in names(ranges.shp.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  validity.ranges.shp.BirdLife[[i]] <- checkValidity(ranges.shp.BirdLife[[i]], quiet=FALSE)
}
rm(i)

## project the polygons into equal-area
for (i in names(ranges.shp.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  for (j in names(ranges.shp.BirdLife[[i]])) {  # loop over taxa
    ranges.shp.BirdLife[[i]][[j]] <- spTransform(ranges.shp.BirdLife[[i]][[j]], CRS=CRS("+init=epsg:3410"))
  }
}
rm(i,j)

## check validity of the equal-area projections (they're identical, so it appears the projection doesn't affect their validity)
validity.ranges.shp.BirdLife.eArea <- list()
for (i in names(ranges.shp.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  validity.ranges.shp.BirdLife.eArea[[i]] <- checkValidity(ranges.shp.BirdLife[[i]], quiet=FALSE)
}
rm(i)

## calculate areas of these equal-area projections
areas.ranges.shp.BirdLife <- list()
for (i in names(ranges.shp.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  areas.ranges.shp.BirdLife[[i]] <- getAreas(ranges.shp.BirdLife[[i]], quiet=FALSE)
}
rm(i)

## make new maps, flattened with gUnaryUnion, to simplify the geometries
ranges.shp.BirdLife.UnaryUnion <- list()
for (i in names(ranges.shp.BirdLife)) {  # loop over year-round, migratory, and/or non-migratory ranges
  ranges.shp.BirdLife.UnaryUnion[[i]] <- list()
  for (j in names(ranges.shp.BirdLife[[i]])) {  # loop over taxa
    ranges.shp.BirdLife.UnaryUnion[[i]][[j]] <- gUnaryUnion(ranges.shp.BirdLife[[i]][[j]])  # generate flattened geometry (all polygon features collapsed into one)
  }
}
rm(i,j)

## calculate areas of these unioned geometries
areas.ranges.shp.BirdLife.UnaryUnion <- list()
for (i in names(ranges.shp.BirdLife.UnaryUnion)) {  # loop over year-round, migratory, and/or non-migratory ranges
    areas.ranges.shp.BirdLife.UnaryUnion[[i]] <- getAreas(ranges.shp.BirdLife.UnaryUnion[[i]], quiet=FALSE)
}
rm(i)

## compare areas of unioned geometries to areas of original equal-area geometries
areas.ranges.shp.BirdLife.changed_by_UnaryUnion <- list()
for (i in names(areas.ranges.shp.BirdLife.UnaryUnion)) {  # loop over year-round, migratory, and/or non-migratory ranges
  areas.ranges.shp.BirdLife.changed_by_UnaryUnion[[i]] <- ((areas.ranges.shp.BirdLife.UnaryUnion[[i]] - areas.ranges.shp.BirdLife[[i]][match(names(areas.ranges.shp.BirdLife.UnaryUnion[[i]]), names(areas.ranges.shp.BirdLife[[i]]))])/areas.ranges.shp.BirdLife[[i]][match(names(areas.ranges.shp.BirdLife.UnaryUnion[[i]]), names(areas.ranges.shp.BirdLife[[i]]))])[which((areas.ranges.shp.BirdLife.UnaryUnion[[i]] - areas.ranges.shp.BirdLife[[i]][match(names(areas.ranges.shp.BirdLife.UnaryUnion[[i]]), names(areas.ranges.shp.BirdLife[[i]]))]) != 0)]
}
rm(i)

areas.ranges.shp.BirdLife.changed_by_UnaryUnion  # many small changes, and some larger ones; all the larger ones are negative


# plot geometries before and after unioning to compare them
for (i in names(areas.ranges.shp.BirdLife.changed_by_UnaryUnion)) {  # loop over year-round, migratory, and non-migratory ranges
  for (j in names(areas.ranges.shp.BirdLife.changed_by_UnaryUnion[[i]])) {  # loop over taxa
    if (abs(areas.ranges.shp.BirdLife.changed_by_UnaryUnion[[i]][[j]]) > 0.0001) {
        cat(i, j, "\n", sep = " ")
        plot(ranges.shp.BirdLife[[i]][[j]], col="red")
        plot(ranges.shp.BirdLife.UnaryUnion[[i]][[j]], col="blue", add=TRUE)
    }
  }
}
rm(i,j)
# check for areas that are red and not covered by blue, as those indicate differences
# checked all such areas in QGIS, and they were due to overlapping features; in all cases where the difference is non-trivial, it's because I lumped two taxa that had overlapping ranges; when flattening those combined ranges, the range size is reduced by the area of the overlap, because it goes from being counted twice (before UnaryUnion) to being counted once (after UnaryUnion)

## check validity of these unioned geometries
validity.ranges.shp.BirdLife.UnaryUnion <- list()
for (i in names(ranges.shp.BirdLife.UnaryUnion)) {  # loop over year-round, migratory, and/or non-migratory ranges
  validity.ranges.shp.BirdLife.UnaryUnion[[i]] <- checkValidity(ranges.shp.BirdLife.UnaryUnion[[i]], quiet=FALSE)
  }
}
rm(i)

## make new maps, buffered by 0, to correct for problems with validity in a few of the geometries; only do this for invalid geometries, to avoid introducing new problems
ranges.shp.BirdLife.UnaryUnion.buffer0 <- ranges.shp.BirdLife.UnaryUnion
for (i in names(ranges.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  for (j in which(!validity.ranges.shp.BirdLife.UnaryUnion[[i]])) {  # loop over taxa with invalid geometries
    print(j)
      ranges.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]] <- gBuffer(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]], width=0)
  }
}
rm(i,j)

## check validity of these unioned, buffered geometries
validity.ranges.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(ranges.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
    validity.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]] <- checkValidity(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]], quiet=FALSE)
}
rm(i)

# print out names of reamining invalid geometries
for (i in names(validity.ranges.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  print(which(!validity.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]]))
}
rm(i)
# no invalid geometries left!

## calculate areas of these unioned, buffered geometries
areas.ranges.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(ranges.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  areas.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]] <- getAreas(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]], quiet=FALSE)
}
rm(i)

## compare areas of unioned geometries to areas of original equal-area geometries
areas.ranges.shp.BirdLife.changed_by_UnaryUnion_buffer0 <- list()
for (i in names(areas.ranges.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  areas.ranges.shp.BirdLife.changed_by_UnaryUnion_buffer0[[i]] <- ((areas.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]] - areas.ranges.shp.BirdLife[[i]][match(names(areas.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(areas.ranges.shp.BirdLife[[i]]))])/areas.ranges.shp.BirdLife[[i]][match(names(areas.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(areas.ranges.shp.BirdLife[[i]]))])[which((areas.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]] - areas.ranges.shp.BirdLife[[i]][match(names(areas.ranges.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(areas.ranges.shp.BirdLife[[i]]))]) != 0)]
}
rm(i)

areas.ranges.shp.BirdLife.changed_by_UnaryUnion_buffer0  # check differences in areas
# it doesn't fill in the holes!!!
# changes are tiny, or are the ones due to lumping overlapping ranges


### calculate overlaps among taxa

picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0 <- list()  # initialize a list to store overlaps
picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[["all_year"]] <- calculateOverlaps(ranges.shp.BirdLife.UnaryUnion.buffer0[["all_year"]], mode="both", quiet=FALSE)  # calculate overlaps treating all ranges as year-round
taxnames.migratory.tmp <- unique(gsub("_summer", "", gsub("_winter", "", names(ranges.shp.BirdLife.UnaryUnion.buffer0[["migratory"]]))))  # generate vector of migratory taxon names so that I can leave them out of the nonmigratory ranges 
picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[["migratory"]] <- calculateOverlapsMigratory(geom.ranges.nonmigratory = ranges.shp.BirdLife.UnaryUnion.buffer0[["all_year"]][!(names(ranges.shp.BirdLife.UnaryUnion.buffer0[["all_year"]]) %in% taxnames.migratory.tmp)], geom.ranges.migratory=ranges.shp.BirdLife.UnaryUnion.buffer0[["migratory"]], mode="both", quiet=FALSE)  # calculate overlaps treating nonmigratory ranges as year-round and migratory ranges split 50% summer / 50% winter
rm(taxnames.migratory.tmp)

## create a version including only the picinae
picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]] <- list()
  for (j in names(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]])) {  # loop over unscaled and scaled overlaps
    picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]] <- picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]][-grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", rownames(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]])), -grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", colnames(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]]))]  # extract overlaps after removing non-Picinae genera
  }
}
rm(i,j)

## calculate the summed overlaps for picidae
picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]] <- list()
  for (j in names(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]])) {  # loop over unscaled and scaled overlaps
    picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]] <- rowSums(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]]) - diag(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]])  # sum the rows, and remove overlap with self (the diagonal of the matrix)
  }
}
rm(i,j)

# calculate the summed overlaps for picinae
picinae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  picinae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]] <- list()
  for (j in names(picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]])) {  # loop over unscaled and scaled overlaps
    picinae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]] <- rowSums(picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]]) - diag(picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]])  # sum the rows, and remove overlap with self (the diagonal of the matrix)
  }
}
rm(i,j)

## calculate range sizes for picidae; need to do this mostly for the migratory species
picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]] <- diag(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][["overlaps"]])  # extract range sizes as the diagonal of the matrix of unscaled overlaps
}
rm(i)

# calculate range sizes for picinae
picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0 <- list()
for (i in names(picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over year-round, migratory, and/or non-migratory ranges
  picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]] <- diag(picinae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[i]][["overlaps"]])  # extract range sizes as the diagonal of the matrix of unscaled overlaps
}
rm(i)

i <- 1
j <- 1
# check that range sizes are the same for nonmigratory taxa and different for migratory taxa
which(!(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]] == picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["all_year"]][match(names(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["all_year"]]))]))  # output range sizes that are different for Picidae in the unscaled overlaps
which(!(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]] == picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["migratory"]][match(names(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["migratory"]]))]))  # output range sizes that are different for Picidae in the scaled overlaps
which(!(picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]] == picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["all_year"]][match(names(picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["all_year"]]))]))  # output range sizes that are different for Picinae in the unscaled overlaps
which(!(picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]] == picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["migratory"]][match(names(picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[i]]), names(picinae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["migratory"]]))]))  # output range sizes that are different for Picinae in the scaled overlaps
rm(i,j)


### save off the range sizes, absolute and scaled overlaps, and the summed absolute and scaled overlaps
save(file='Picidae_overlaps.RData', list=c(grep("rangesize.shp", ls(), value=TRUE), grep("overlaps.shp", ls(), value=TRUE)))


### plot histograms of overlaps
hist(log(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["all_year"]]))
hist(log(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[["migratory"]]))

hist(log(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[["all_year"]][["overlaps"]]))
hist(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[["all_year"]][["overlaps.scaled"]])


### mapping ranges for species richness map (Appendix 4)
## I treat ranges of migratory species as being occupied year-round, as weighting them is difficult for the purposes of this map

i <- "all_year"

# get the limits of the ranges so that I can build an appropriate bounding box
picidae.bbox.by_species <- array(dim=c(2,2,length(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]])))
for (j in 1:length(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]])) {
  picidae.bbox.by_species[,,j] <- bbox(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]])
}

# generate a universal bounding box using the minima and maxima from the species-level bounding boxes
picidae.bbox <- matrix(nrow=2, ncol=2, dimnames=list(c("x","y"), c("min","max")))
picidae.bbox[1,1] <- min(picidae.bbox.by_species[1,1,])
picidae.bbox[2,1] <- min(picidae.bbox.by_species[2,1,])
picidae.bbox[1,2] <- max(picidae.bbox.by_species[1,2,])
picidae.bbox[2,2] <- max(picidae.bbox.by_species[2,2,])

picidae.bbox.spatial <- matrix(nrow=5, ncol=2)
picidae.bbox.spatial[1,] <- picidae.bbox[,1]
picidae.bbox.spatial[2,] <- c(picidae.bbox[1,2], picidae.bbox[2,1])
picidae.bbox.spatial[3,] <- picidae.bbox[,2]
picidae.bbox.spatial[4,] <- c(picidae.bbox[1,1], picidae.bbox[2,2])
picidae.bbox.spatial[5,] <- picidae.bbox[,1]

# make the bounding box a SpatialPolygons object
picidae.bbox.spatial <- Polygon(picidae.bbox.spatial)
picidae.bbox.spatial <- Polygons(list(picidae.bbox.spatial),1)
picidae.bbox.spatial <- SpatialPolygons(list(picidae.bbox.spatial))
proj4string(picidae.bbox.spatial) = CRS("+init=epsg:3410")

# start a rangeMapper session
picidae.rangeMap <- rangeMap.start(file="picidae_rangeMap.sqlite", dir="/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Picidae_distribution_data/", overwrite=TRUE)
global.bbox.save(con=picidae.rangeMap, bbox=picidae.bbox.spatial)


ranges.shp.BirdLife.UnaryUnion.buffer0.polygons <- list()
for (j in 1:length(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]])) {
  ranges.shp.BirdLife.UnaryUnion.buffer0.polygons[[j]] <- ranges.shp.BirdLife.UnaryUnion.buffer0[[i]][[j]]@polygons
  ranges.shp.BirdLife.UnaryUnion.buffer0.polygons[[j]][[1]]@ID <- names(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]])[j]
}

ranges.shp.BirdLife.UnaryUnion.buffer0.SpatPoly <- SpatialPolygons(Srl=sapply(ranges.shp.BirdLife.UnaryUnion.buffer0.polygons, function(x) x[[1]]), pO=1:length(ranges.shp.BirdLife.UnaryUnion.buffer0.polygons), proj4string=CRS("+init=epsg:3410"))

ranges.shp.BirdLife.UnaryUnion.buffer0.DF.placeholder <- data.frame("taxon"=names(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]]))
rownames(ranges.shp.BirdLife.UnaryUnion.buffer0.DF.placeholder) <- names(ranges.shp.BirdLife.UnaryUnion.buffer0[[i]])

ranges.shp.BirdLife.UnaryUnion.buffer0.SpatPolyDF <- SpatialPolygonsDataFrame(Sr=ranges.shp.BirdLife.UnaryUnion.buffer0.SpatPoly, data=ranges.shp.BirdLife.UnaryUnion.buffer0.DF.placeholder)

## save a grid size and canvas to the rangeMapper object
gridSize.save(con=picidae.rangeMap, gridSize=(picidae.bbox[2,2]-picidae.bbox[2,1])/500)
canvas.save(con=picidae.rangeMap)


processRanges(con=picidae.rangeMap, spdf=ranges.shp.BirdLife.UnaryUnion.buffer0.SpatPolyDF, ID="taxon")  # quantify overlap in ranges
rangeMap.save(picidae.rangeMap)  # save rangeMap

picidae.rangeMap.species_richness <- rangeMap.fetch(picidae.rangeMap)  # retrieve the rangeMap so that it doesn't need to be recreated in later instances

## output the map of species richness
pdf(file="species_richness_map.pdf", width=12, height=6)
plot(picidae.rangeMap.species_richness)
dev.off()

## output the map with grayscale
pdf(file="species_richness_map_2.pdf", width=12, height=6)
spplot(rangeMap.fetch(picidae.rangeMap), col.regions=gray.colors(n=20))
dev.off()

## import boundaries of landforms (continents, islands, large lakes, etc.)
landforms <- readOGR(dsn="ne_110m_land/", layer="ne_110m_land")
landforms <- spTransform(landforms, CRS("+init=epsg:3410"))

## output map with landforms shown
pdf(file="species_richness_map_3.pdf", width=9.5, height=6)
spplot(picidae.rangeMap.species_richness, col.regions=rev(gray.colors(n=16, start=0, end=0.9, gamma=1)), panel=function(...) {
  panel.levelplot(...)
  sp.polygons(landforms, fill=NA, col="black", lwd=0.4)
})
dev.off()
