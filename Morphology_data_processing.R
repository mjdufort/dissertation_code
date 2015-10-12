## This is is one of several files containing scripts and functions used in processing and analysis of data for Matthew Dufort's Ph.D. dissertation at the University of Minnesota, titled "Coexistence, Ecomorphology, and Diversification in the Avian Family Picidae (Woodpeckers and Allies)."

## this file contains functions and scripts to process the morphological data to the point that it can be used in my community-level (Chapter 2) and broader comparative (Chapter 3) analyses

## load necessary packages
library(plyr)
library(nlme)
library(ape)
library(geiger)
library(phytools)


### load trees, overlaps, and other data, and generate basic functions

## load trees from Tree_manipulation.R
load(file='R_trees_for_combined_analyses.RData')
picidae.RAxML.all.BEAST_calibrated[grep("rate", names(picidae.RAxML.all.BEAST_calibrated))] <- NULL  # remove extra stuff that I don't need from the BEAST-calibrated RAxML tree
picidae.RAxML.all.BEAST_calibrated.no_calibs[grep("rate", names(picidae.RAxML.all.BEAST_calibrated.no_calibs))] <- NULL  # remove extra stuff that I don't need from the BEAST-calibrated RAxML tree

## drop outgroups from trees
taxa.remove.final_trees <- c("Indicatoridae_spp", "Ramphastidae_spp")
picidae.RAxML.all.BEAST_calibrated <- drop.tip(picidae.RAxML.all.BEAST_calibrated, tip=taxa.remove.final_trees)
picidae.RAxML.all.BEAST_calibrated.no_calibs <- drop.tip(picidae.RAxML.all.BEAST_calibrated.no_calibs, tip=taxa.remove.final_trees)
picidae.RAxML.all.r8s_calibrated.ladderized <- drop.tip(picidae.RAxML.all.r8s_calibrated.ladderized, tip=taxa.remove.final_trees)
picidae.RAxML.all.r8s_calibrated.no_calibs.ladderized <- drop.tip(picidae.RAxML.all.r8s_calibrated.no_calibs.ladderized, tip=taxa.remove.final_trees)
picidae.RAxML.all.rooted.ladderized <- drop.tip(picidae.RAxML.all.rooted.ladderized, tip=taxa.remove.final_trees)
picidae.RAxML.all.rooted.ladderized.chronos.lambda1 <- drop.tip(picidae.RAxML.all.rooted.ladderized.chronos.lambda1, tip=taxa.remove.final_trees)
picidae.mrbayes.all.ladderized <- drop.tip(picidae.mrbayes.all.ladderized, tip=taxa.remove.final_trees)
rm(taxa.remove.final_trees)

## load the absolute and scaled overlaps from Quantifying_distribution_overlaps.R
load(file='Picidae_overlaps.RData')

## the function reduceToMeanByCols() is a wrapper for ddply to reduce a data frame by a set of factor or character variables (reduce_colnames), taking the means of another set of variables (column numbers in mean_cols)
# as input, it takes data (a data frame), reduce_colnames (the names of the columns to reduce by), and mean_cols (the columns to take means of over the factors in reduce_colnames)
# it returns a data frame reduced to having only one row for each combination of the values in reduce_colnames
reduceToMeanByCols <- function(data, reduce_colnames, mean_cols) {
  library(plyr)
  data.reduced <- ddply(.data=data, .variables=reduce_colnames, function(x){
    y <- subset(x, select=mean_cols)
    apply(y, 2, mean, na.rm=TRUE)
  })
  return(data.reduced)
}


###### run the analyses for Picidae ######

### read in the morphological data from the database (caliper measurements)
picidae.morph <- list()
picidae.morph$raw <- read.csv('Picidae_morphology_data.csv', header=T, stringsAsFactors=FALSE)  # read in the raw morphology data, exported from my database
picidae.morph$raw$Data_date <- as.Date(picidae.morph$raw$Data_date, format="%m/%d/%y")  # convert the date strings to date objects

picidae.inst_catnum.exclude <- c("KU_96827","USNM_431833","USNM_431677","USNM_347456","BellMNH_11400","USNM_621108","UMMZ_74129","USNM_430385","USNM_431832","FLMNH_27223","FLMNH_47211","BellMNH_42345","UMMZ_224192","UWBM_42166") # exclude some specimens with questionable data

picidae.morph$raw.exclude_bad <- picidae.morph$raw[!((picidae.morph$raw$My_genus_species=="") | (picidae.morph$raw$Inst_CatNum %in% picidae.inst_catnum.exclude)),]  # exclude rows with blank My_genus_species and those I selected to exclude
picidae.morph$raw.exclude_bad_juvs <- picidae.morph$raw.exclude_bad[-grep(pattern="juv", x=picidae.morph$raw.exclude_bad$AgeClass, ignore.case=TRUE),]  # exclude those with "juv" in the AgeClass field

picidae.morph$log <- cbind(picidae.morph$raw.exclude_bad_juvs[1:14],log(picidae.morph$raw.exclude_bad_juvs[15:ncol(picidae.morph$raw.exclude_bad_juvs)]))  # log-transform all measurements


### check data for problems and outliers

## plot histograms of all measurements for each species
picidae.species.names <- sort(unique(picidae.morph$log$My_genus_species))  # extract the taxon names as a character vector
for (i in picidae.species.names) {  # loop over taxon
  pdf(file=paste(gsub("/", "_", i), "_histogram_matrix.pdf", sep=""), width=10, height=8)
  par(mfrow=c(4,5))
  for (j in 15:34) {  # loop over measurements
    data.tmp <- picidae.morph$log[picidae.morph$log$My_genus_species==i, j]  # pull out the data for the current species
    if (any(!is.na(data.tmp))) {hist(data.tmp, main=colnames(picidae.morph$log)[j], breaks=10)  # if there is data for that measurement, plot a histogram
    } else plot.new()
  }
  dev.off()
  rm(data.tmp)
}
rm(i,j)

## generate histogram of Melanerpes carolinus without the repeated measurements of my reference specimen
pdf(file="Melanerpes_carolinus_ex_reference_histogram_matrix.pdf", width=10, height=8)
par(mfrow=c(4,5))
for (j in 15:34) {
  hist(picidae.morph$log[((picidae.morph$log$My_genus_species=="Melanerpes_carolinus") & (picidae.morph$log$Inst_CatNum != "BellMNH_47191")), j], main=colnames(picidae.morph$log)[j], breaks=10)
  }
dev.off()

## generate histogram of Sphyrapicus varius without the repeated measurements of my reference specimen
pdf(file="Sphyrapicus_varius_ex_reference_histogram_matrix.pdf", width=10, height=8)
par(mfrow=c(4,5))
for (j in 15:34) {
  data.tmp <- picidae.morph$log[((picidae.morph$log$My_genus_species=="Sphyrapicus_varius") & (picidae.morph$log$Inst_CatNum != "BellMNH_47202")), j]
  if (any(!is.na(data.tmp))) {hist(data.tmp, main=colnames(picidae.morph$log)[j], breaks=10)
  } else plot.new()
}
dev.off()
rm(data.tmp)

## generate scatterplot matrix for measurements of Colaptes punctigula, to visualize weird differences in size/shape among different measurements
pdf(file="/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Picidae_combined_analyses_dist_morph_evol_divers/Picidae_morph_data_checking/Colaptes_punctigula_scatterplot_matrix.pdf", width=40, height=40)
par(mfrow=c(20,20))
data.tmp <- picidae.morph$log[picidae.morph$log$My_genus_species=="Colaptes_punctigula",15:34]  # extract measurement data for Colaptes punctigula
for (i in 1:20) {
  for (j in 1:20) {
    if (any (!(is.na(data.tmp[,i]) | is.na(data.tmp[,j])))) {
      plot(data.tmp[,i], data.tmp[,j], xlim=c(min(data.tmp[,i]-(0.1*(range(data.tmp[,i], na.rm=TRUE)[2] - range(data.tmp[,i], na.rm=TRUE)[1])), na.rm=TRUE), max(data.tmp[,i]+(0.1*(range(data.tmp[,i], na.rm=TRUE)[2] - range(data.tmp[,i], na.rm=TRUE)[1])), na.rm=TRUE)), ylim=c(min(data.tmp[,j]-(0.1*(range(data.tmp[,j], na.rm=TRUE)[2] - range(data.tmp[,j], na.rm=TRUE)[1])), na.rm=TRUE), max(data.tmp[,j]+(0.1*(range(data.tmp[,j], na.rm=TRUE)[2] - range(data.tmp[,j], na.rm=TRUE)[1])), na.rm=TRUE)), xlab=colnames(data.tmp)[i], ylab=colnames(data.tmp)[j])
    } else plot.new()
  }
}
dev.off()
rm(data.tmp)


### reduce the data by left/right and extra variables
## generate columns with means for data that have left and right
picidae.morph$log$Mandible_mean_ramp_present <- apply(cbind(picidae.morph$log$Mandible_left_ramp_present,picidae.morph$log$Mandible_right_ramp_present),MARGIN=1,mean,na.rm=TRUE)
picidae.morph$log$Mandible_mean_ramp_absent <- apply(cbind(picidae.morph$log$Mandible_left_ramp_absent,picidae.morph$log$Mandible_right_ramp_absent),MARGIN=1,mean,na.rm=TRUE)
picidae.morph$log$Humerus_mean <- apply(cbind(picidae.morph$log$Humerus_left,picidae.morph$log$Humerus_right),MARGIN=1,mean,na.rm=TRUE)
picidae.morph$log$Ulna_mean <- apply(cbind(picidae.morph$log$Ulna_left,picidae.morph$log$Ulna_right),MARGIN=1,mean,na.rm=TRUE)
picidae.morph$log$Femur_mean <- apply(cbind(picidae.morph$log$Femur_left,picidae.morph$log$Femur_right),MARGIN=1,mean,na.rm=TRUE)
picidae.morph$log$Tibiotarsus_mean <- apply(cbind(picidae.morph$log$Tibiotarsus_left,picidae.morph$log$Tibiotarsus_right),MARGIN=1,mean,na.rm=TRUE)
picidae.morph$log$Tarsometatarsus_mean <- apply(cbind(picidae.morph$log$Tarsometatarsus_left,picidae.morph$log$Tarsometatarsus_right),MARGIN=1,mean,na.rm=TRUE)

picidae.morph.log.reduced_var <- list()  # intiialize a list to store morphological data with reduced number of variables
picidae.morph.log.reduced_var$base <- subset(picidae.morph$log, select=-c(1,7:14, grep("left|right",names(picidae.morph$log)))) # strip out the columns from the morphology data that I'm going to use; this removes the extraneous non-measurement variables, and the variables with left or right in the name
picidae.morph.log.reduced_var$base <- cbind(picidae.morph.log.reduced_var$base[-grep("Pygo", names(picidae.morph.log.reduced_var$base))], picidae.morph.log.reduced_var$base[grep("Pygo", names(picidae.morph.log.reduced_var$base))])  # re-order the columns so that pygostyle is at the end


### read in the raw data from the measurements from images

picidae.morph_from_images <- list()  # initialize a list to store measurements
picidae.morph_from_images$raw <- read.csv('Picidae_ImageJ_measurements.csv', header=T, stringsAsFactors=FALSE) # read in the raw morphology data taken from images, exported from my Excel file
picidae.morph_from_images$log <- cbind(picidae.morph_from_images$raw[c(1:3,8:9)], log(picidae.morph_from_images$raw[,4:7])) # log-transform the measuremetn values
picidae.morph_from_images$log <- picidae.morph_from_images$log[!(picidae.morph_from_images$log$Inst_CatNum %in% c(picidae.inst_catnum.exclude, picidae.morph$raw[grep(pattern="juv", x=picidae.morph$raw$AgeClass, ignore.case=TRUE),"Inst_CatNum"], picidae.morph$raw[(picidae.morph$raw$My_genus_species==""),"Inst_CatNum"])),] # exclude individuals with bad data, and juveniles, and those with a blank My_genus_species

## generate histograms of the maxilla measurements from images
for (i in picidae.species.names) {  # loop over taxon
  pdf(file=paste(gsub("/", "_", i), "_image_histogram_matrix.pdf", sep=""), width=5, height=5)
  par(mfrow=c(2,2))
  for (j in 6:9) {  # loop over measurements
    data.tmp <- picidae.morph_from_images$log[picidae.morph_from_images$log$Inst_CatNum %in% picidae.morph$log$Inst_CatNum[picidae.morph$log$My_genus_species == i], j]  # extract data for current taxon and measurement
    if (any(!is.na(data.tmp))) {hist(data.tmp, main=colnames(picidae.morph_from_images$log)[j], breaks=10)  # plot histogram if there is any data present
    } else plot.new()
  }
  dev.off()
  rm(data.tmp)
}
rm(i,j)

picidae.morph_from_images$log <- picidae.morph_from_images$log[-grep("depth", colnames(picidae.morph_from_images$log))]  # strip out maxilla_depth, because it's not very reproducible
picidae.morph_from_images$log.reduced_inds <- reduceToMeanByCols(picidae.morph_from_images$log, reduce_colnames="Inst_CatNum", mean_cols=6:ncol(picidae.morph_from_images$log))  # calculate averages by individual where there are multiple measurements of a single individual


### quantify date variation in calipers and adjust for it


## the function caliperVariationSumSq() calculates the summed squared deviations from expected values of measurements, given a slope, intercept, and a date-by-date adjustment factor for the slope or intercept; it is used for optimizing the values of these parameters
# as input, it takes pars (a vector containing the proposed values of slope or intercept, and the values of date.adjust, which are date-specific intercepts or slopes), obs (a matrix of measurement values for a single individual, with dates), and fitmethod (if fitmethod="slope", it assumes pars[1] is a universal intercept and pars[2:length(pars)] are date-specific slopes; if fitmethod="intercept", it assumes pars[1] is a universal slope and pars[2:length(pars)] are date-specific intercepts)
# both the data in obs and the parameter values in pars should be sorted in the same order
# in either case, the universal parameter is pars[1], and the date-specific parameters are pars[2:length(pars)]
caliperVariationSumSq <- function(pars, obs, fitmethod, quiet=TRUE) {
  obs.mean <- matrix(apply(obs, MARGIN=2, FUN=mean, na.rm=TRUE), nrow=1)  # calculate the mean values of each variable
  date.adjust <- as.matrix(pars[2:length(pars)], byrow=TRUE)  # extract the date-specific parameter values
  if (fitmethod == "slope") {
    preds <- (date.adjust %*% obs.mean) + pars[1]  # calculate predicted values with universal intercept and date-specific slopes
  }
  if (fitmethod == "intercept") {
    preds <- matrix(rep((pars[1] * obs.mean), times=nrow(obs)), ncol=ncol(obs), byrow=TRUE) + matrix(rep(date.adjust, times=ncol(obs)), nrow=nrow(obs), byrow=FALSE)  # calculate predicted values with universal slope and date-specific intercepts
  }
  
  sumsq <- sum(((preds - obs)^2), na.rm=TRUE)  # calculate and sum the squared deviations
  if (!quiet) print(sumsq)
  return(sumsq)
}

## the function optimizeCaliperVariation() finds the optimum values of a universal slope or intercept and date-specific intercepts or slopes
# as input, it takes obs.with_dates (a matrix of measurements from a single individual with the dates of measurements in the first column), and fitmethod (if fitmethod="slope", it assumes pars[1] is a universal intercept and pars[2:length(pars)] are date-specific slopes; if fitmethod="intercept", it assumes pars[1] is a universal slope and pars[2:length(pars)] are date-specific intercepts)
# it returns a list containing the method used ("slope" or "intercept"), the value of the fixed parameter, the values of the date-specific parameter, and the sum of squared deviations using the optimized parameter values
optimizeCaliperVariation <- function(obs.with_dates, fitmethod, quiet=TRUE) {
  dates <- obs.with_dates[,1]  # extract the dates from the data matrix
  obs.no_dates <- obs.with_dates[,-1]  # extract data without dates from the data matrix
  
  ## set initial values of parameters
  if (fitmethod=="intercept") {inits <- c(1,rep(0,times=nrow(obs.no_dates)))}
  if (fitmethod=="slope") {inits <- c(0, rep(1,times=nrow(obs.no_dates)))}
  
  date.adjustments.optimized <- optim(inits, fn=caliperVariationSumSq, obs=obs.no_dates, fitmethod=fitmethod, quiet=quiet, method="BFGS")  # optimize the parameters
  if (!quiet) print(date.adjustments.optimized)
  date.adjusts <- date.adjustments.optimized$par[2:length(date.adjustments.optimized$par)]  # extract the date-specific parameter values from the list returned by parameter optimization
  names(date.adjusts) <- dates  # assign dates as names for the elements of the vector containing date-specific parameters
  return(list(fitmethod=fitmethod, param.fixed=date.adjustments.optimized$par[1], date.adjusts=date.adjusts, sumsq=date.adjustments.optimized$value))
}

obs.Mel_car.reference <- picidae.morph.log.reduced_var$base[picidae.morph.log.reduced_var$base$Inst_CatNum=="BellMNH_47191",-(2:5)] # subset the data for the columns I want and only the reference Melanerpes carolinus
caliper.variation.optimized_slope.Mel_car.reference <- optimizeCaliperVariation(obs.Mel_car.reference, fitmethod="slope", quiet=FALSE)  # optimize parameter values with fixed intercept and date-specific slope
caliper.variation.optimized_intercept.Mel_car.reference <- optimizeCaliperVariation(obs.Mel_car.reference, fitmethod="intercept", quiet=FALSE)  # optimize parameter values with fixed slope and date-specific intercept

obs.Sph_var.reference <- picidae.morph.log.reduced_var$base[picidae.morph.log.reduced_var$base$Inst_CatNum=="BellMNH_47202",-(2:5)] # subset the data for the columns I want and only the reference Sphyrapicus varius
caliper.variation.optimized_slope.Sph_var.reference <- optimizeCaliperVariation(obs.Sph_var.reference, fitmethod="slope", quiet=FALSE)  # optimize parameter values with fixed intercept and date-specific slope
caliper.variation.optimized_intercept.Sph_var.reference <- optimizeCaliperVariation(obs.Sph_var.reference, fitmethod="intercept", quiet=FALSE)  # optimize parameter values with fixed slope and date-specific intercept

## for both reference specimens, sumsq is definitely lower with date-specific slopes!; date adjust values are fairly similar between Mel_car and Sph_var when the adjustment values are large (which is when they matter)

## calculate average intercept and the average slopes by date from Sph_var and Mel_car reference measurements
caliper.variation.optimized_slope.averaged <- list(intercept=mean(c(caliper.variation.optimized_slope.Sph_var.reference$param.fixed, caliper.variation.optimized_slope.Mel_car.reference$param.fixed)), date.slopes=numeric(length(unique(c(names(caliper.variation.optimized_slope.Sph_var.reference$date.adjusts), names(caliper.variation.optimized_slope.Mel_car.reference$date.adjusts))))))  # create a list containing the average of the fixed parameter value and a vector to store the date-specific parameter values
names(caliper.variation.optimized_slope.averaged$date.slopes) <- as.character(sort(unique(c(names(caliper.variation.optimized_slope.Sph_var.reference$date.adjusts), names(caliper.variation.optimized_slope.Mel_car.reference$date.adjusts)))))  # assign dates with reference measurements as names to the vector of date-specific parameter values

for (i in 1:length(caliper.variation.optimized_slope.averaged$date.slopes)) {caliper.variation.optimized_slope.averaged$date.slopes[i] <- mean(c(caliper.variation.optimized_slope.Sph_var.reference$date.adjusts[names(caliper.variation.optimized_slope.Sph_var.reference$date.adjusts)==names(caliper.variation.optimized_slope.averaged$date.slopes)[i]], caliper.variation.optimized_slope.Mel_car.reference$date.adjusts[names(caliper.variation.optimized_slope.Mel_car.reference$date.adjusts)==names(caliper.variation.optimized_slope.averaged$date.slopes)[i]]), na.rm=TRUE)}  # loop over the dates with reference measurements, calculating the arithmetic mean value of the date-specific parameter
rm(i)

rm(obs.Mel_car.reference, obs.Sph_var.reference, caliper.variation.optimized_intercept.Mel_car.reference, caliper.variation.optimized_intercept.Sph_var.reference, caliper.variation.optimized_slope.Mel_car.reference, caliper.variation.optimized_slope.Sph_var.reference)  # remove unneeded variables from this section

### next step is to  generate adjustments for the data, and apply them to the data set; this requires calculating some second-order adjustments, as I don't have reference measurements for every date
## the set of functions below uses deviations between measurements of the same specimen on different takes to find optimal adjustments to intercept or slope for a date without adjustments based on reference specimens

## the function interceptSumsq() calculates the sum of squared deviations for use in optimizing the intercept for date adjustment of data
# as input, it takes intercept.new (the proposed intercept for the new date), intercept.orig (the intercept for the date with adjustment already determined), slope (the slope), data (the data as a data frame, including only specimens for which measurements are present for both dates), date.new (the date for which we're calculating a new intercept), and date.orig (the date for which we have an intercept)
# it returns the sum of the squared deviations between the data for each date, for specimens with measurements for both dates, after adjusting data from both dates using the given values of intercepts and slope. 
interceptSumsq <- function(intercept.new, intercept.orig, slope, data, date.new, date.orig) {
  sums <- 0  # intialize a numeric to store the sum or squared differences
  Inst_CatNums <- unique(data$Inst_CatNum)  # extract a vector of individual IDs
  for (i in 1:length(Inst_CatNums)) {  # loop over individuals
    sum.tmp <- sum((((as.matrix(data[data$Inst_CatNum==Inst_CatNums[i] & data$Data_date==date.orig, -c(1:2)]) - intercept.orig) / slope) - (as.matrix((data[data$Inst_CatNum==Inst_CatNums[i] & data$Data_date==date.new, -c(1:2)]) - intercept.new) / slope))^2, na.rm=TRUE)  # for the current individual, adjust the data for each date, take the difference, square the difference, and sum over all measurements
    sums <- sums + sum.tmp  # add the sum for the current individual to the overall sum
  }
  return(sums)
}

## the function slopeSumsq() calculates the sum of squared deviations for use in optimizing the slope for date adjustment of data
# inputs and value returned are same as for interceptSumsq(), except that there is a single intercept and date-specific slopes
slopeSumsq <- function(slope.new, slope.orig, intercept, data, date.new, date.orig) {
  sums <- 0  # intialize a numeric to store the sum or squared differences
  Inst_CatNums <- unique(data$Inst_CatNum)  # extract a vector of individual IDs
  for (i in 1:length(Inst_CatNums)) {  # loop over individuals
    sum.tmp <- sum((((as.matrix(data[data$Inst_CatNum==Inst_CatNums[i] & data$Data_date==date.orig, -c(1:2)]) - intercept) / slope.orig) - (as.matrix((data[data$Inst_CatNum==Inst_CatNums[i] & data$Data_date==date.new, -c(1:2)]) - intercept) / slope.new))^2, na.rm=TRUE)  # for the current individual, adjust the data for each date, take the difference, square the difference, and sum over all measurements
    sums <- sums + sum.tmp  # add the sum for the current individual to the overall
  }
  return(sums)
}

## the function interceptNewdate() calculates the best intercept for a new date using data from specimens measured on both the new date and a date for which a calculated intercept already exists
# as input, it takes intercept.orig (the calculated intercept for the date with adjustment already determined), slope (the value of the slope), data (the data as a data frame, including only specimens for which measurements are present for both dates), date.new (the date for which an intercept is to be calculated), and date.orig (the date for which an intercept has already been calculated)
# it returns a list containing the optimized value of the intercept for use in adjustment of data for date variation, and the sum of squared deviations using the optimum value
interceptNewdate <- function(intercept.orig, slope, data, date.new, date.orig) {
  optimized.intercept.new <- optimize(f=interceptSumsq, interval=c(-0.1,0.1), intercept.orig=intercept.orig, slope=slope, data=data, date.new=date.new, date.orig=date.orig)  # use the interceptSumsq function to optimize the intercept for the new date
  return(list(intercept.optimized=optimized.intercept.new$minimum, sumsq=optimized.intercept.new$objective))
}

## the function slopeNewdate() calculates the best slope for a new date using data from specimens measured on both the new date and a date for which a calculated slope already exists
# inputs and value returned are same as for interceptNewdate(), except that there is a single intercept and date-specific slopes
slopeNewdate <- function(slope.orig, intercept, data, date.new, date.orig) {
  optimized.slope.new <- optimize(f=slopeSumsq, interval=c(-0.99,1.01), slope.orig=slope.orig, intercept=intercept, data=data, date.new=date.new, date.orig=date.orig)
  return(list(slope.optimized=optimized.slope.new$minimum, sumsq=optimized.slope.new$objective))  # use the slopeSumsq function to optimize the slope for the new date
}


## now calculate the adjustments
caliper.adjust <- list(intercept=caliper.variation.optimized_slope.averaged$intercept, date.slopes=numeric(length(unique(picidae.morph.log.reduced_var$base$Data_date)))) # generate an empty object to store the slope and the date-specific intercepts for ALL dates
names(caliper.adjust$date.slopes) <- sort(unique(picidae.morph.log.reduced_var$base$Data_date)) # set the names for the date-specific slopes
caliper.adjust$date.slopes[names(caliper.variation.optimized_slope.averaged$date.slopes)] <- caliper.variation.optimized_slope.averaged$date.slopes # set the first-order adjustments (where I have reference measurements from that date)

## estimate and fill in the date-specific slopes for dates with/without reference and with overlapping specimens; do this iteratively, until I no longer have any such overlaps

# round 1; using first-order intersections of data with dates reference data
dates.no_reference.1 <- names(caliper.adjust$date.slopes[caliper.adjust$date.slopes==0])  # get dates without calculated date-specific slopes
dates.with_reference.1 <- names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])  # get dates with calculated date-specific slopes
inst_catnum.with_reference.1 <- unique(picidae.morph.log.reduced_var$base$Inst_CatNum[as.character(picidae.morph.log.reduced_var$base$Data_date) %in% dates.with_reference.1])  # get specimens measured on dates with calculated date-specific slopes
inst_catnum.no_reference.1 <- unique(picidae.morph.log.reduced_var$base$Inst_CatNum[as.character(picidae.morph.log.reduced_var$base$Data_date) %in% dates.no_reference.1])  # get specimens measured on dates without calculated date-specific slopes
inst_catnum.overlapping.1 <- inst_catnum.no_reference.1[inst_catnum.no_reference.1 %in% inst_catnum.with_reference.1]  # identify specimens measured on dates with and without calculated date-specific slopes
picidae.morph.log.reduced_var$base[picidae.morph.log.reduced_var$base$Inst_CatNum %in% inst_catnum.overlapping.1,]  # output data for those specimens; use the dates from this output to identify dates for the next step

# first date
slope.optimized.date1 <- slopeNewdate(slope.orig=caliper.adjust$date.slopes["2013-10-21"], intercept=caliper.adjust$intercept, data=picidae.morph.log.reduced_var$base[(picidae.morph.log.reduced_var$base$Inst_CatNum %in% inst_catnum.overlapping.1) & picidae.morph.log.reduced_var$base$Data_date %in% as.Date(c("2011-08-04", "2013-10-21")), -c(3:5)], date.new=as.Date("2011-08-04"), date.orig=as.Date("2013-10-21"))  # find the optimal slope for the new date, as the slope that minimizes the squared deviations between the adjusted values for the measurements from the two dates
caliper.adjust$date.slopes["2011-08-04"] <- slope.optimized.date1$slope.optimized  # store the calculated slope in my vector containing all date-specific slopes

# second date
slope.optimized.date2 <- slopeNewdate(slope.orig=caliper.adjust$date.slope["2013-09-11"], intercept=caliper.adjust$intercept, data=picidae.morph.log.reduced_var$base[(picidae.morph.log.reduced_var$base$Inst_CatNum %in% inst_catnum.overlapping.1) & picidae.morph.log.reduced_var$base$Data_date %in% as.Date(c("2012-06-07", "2013-09-11")), -c(3:5)], date.new=as.Date("2012-06-07"), date.orig=as.Date("2013-09-11"))  # find the optimal slope for the new date, as the slope that minimizes the squared deviations between the adjusted values for the measurements from the two dates
caliper.adjust$date.slopes["2012-06-07"] <- slope.optimized.date2$slope.optimized  # store the calculated slope in my vector containing all date-specific slopes

# round 2; using second-order intersections of data with dates reference data
dates.no_reference.2 <- names(caliper.adjust$date.slopes[caliper.adjust$date.slopes==0])  # get dates without calculated date-specific slopes
dates.with_reference.2 <- names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])  # get dates with calculated date-specific slopes
inst_catnum.with_reference.2 <- unique(picidae.morph.log.reduced_var$base$Inst_CatNum[as.character(picidae.morph.log.reduced_var$base$Data_date) %in% dates.with_reference.2])  # get specimens measured on dates with calculated date-specific slopes
inst_catnum.no_reference.2 <- unique(picidae.morph.log.reduced_var$base$Inst_CatNum[as.character(picidae.morph.log.reduced_var$base$Data_date) %in% dates.no_reference.2])  # get specimens measured on dates without calculated date-specific slopes
inst_catnum.overlapping.2 <- inst_catnum.no_reference.2[inst_catnum.no_reference.2 %in% inst_catnum.with_reference.2]  # identify specimens measured on dates with and without calculated date-specific slopes
picidae.morph.log.reduced_var$base[picidae.morph.log.reduced_var$base$Inst_CatNum %in% inst_catnum.overlapping.2,]  # output data for those specimens; use the dates from this output to identify dates for the next step

# third date
slope.optimized.date3 <- slopeNewdate(slope.orig=caliper.adjust$date.slopes["2012-06-07"], intercept=caliper.variation.optimized_slope.averaged$intercept, data=picidae.morph.log.reduced_var$base[(picidae.morph.log.reduced_var$base$Inst_CatNum %in% inst_catnum.overlapping.2) & picidae.morph.log.reduced_var$base$Data_date %in% as.Date(c("2011-12-09","2012-06-07")), -c(3:5,10:11)], date.new=as.Date("2011-12-09"), date.orig=as.Date("2012-06-07"))  # find the optimal slope for the new date, as the slope that minimizes the squared deviations between the adjusted values for the measurements from the two dates
caliper.adjust$date.slopes["2011-12-09"] <- slope.optimized.date3$slope.optimized  # store the calculated slope in my vector containing all date-specific slopes

# round 3; using third-order intersections of data with dates reference data
dates.no_reference.3 <- names(caliper.adjust$date.slopes[caliper.adjust$date.slopes==0])  # get dates without calculated date-specific slopes
dates.with_reference.3 <- names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])  # get dates with calculated date-specific slopes
inst_catnum.with_reference.3 <- unique(picidae.morph.log.reduced_var$base$Inst_CatNum[as.character(picidae.morph.log.reduced_var$base$Data_date) %in% dates.with_reference.3])  # get specimens measured on dates with calculated date-specific slopes
inst_catnum.no_reference.3 <- unique(picidae.morph.log.reduced_var$base$Inst_CatNum[as.character(picidae.morph.log.reduced_var$base$Data_date) %in% dates.no_reference.3])  # get specimens measured on dates without calculated date-specific slopes
inst_catnum.overlapping.3 <- inst_catnum.no_reference.3[inst_catnum.no_reference.3 %in% inst_catnum.with_reference.3]  # identify specimens measured on dates with and without calculated date-specific slopes
picidae.morph.log.reduced_var$base[picidae.morph.log.reduced_var$base$Inst_CatNum %in% inst_catnum.overlapping.3,]  # output data for those specimens; use the dates from this output to identify dates for the next step

# no more matching dates

## finally, estimate the date-specific slope for dates that don't have overlapping specimens with dates with reference data
# I did this by fitting a loess curve, with the span determined by replotting the model and checking it visually, then using the predicted values from that loess fit for the dates with no intercept values
slope.by.date.loess <- loess(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0] ~ as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0]))), span=0.45, family="gaussian")  # fit a loess curve to the date-specific slopes, using the date as the indepdendent variable and the slope as the dependent variable
plot(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0] ~ as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0]))))  # plot the points with already calculated date-specific slopes
lines(min(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])))):max(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])))), predict(slope.by.date.loess, min(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])))):max(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0]))))))  # plot the loess curve

caliper.adjust$date.slopes[caliper.adjust$date.slopes==0] <- predict(slope.by.date.loess, newdata=as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes==0]))))  # replace empty values with the predicted values from the loess fit
caliper.adjust$date.slopes[is.na(caliper.adjust$date.slopes)] <- caliper.adjust$date.slopes[min(which(!is.na(caliper.adjust$date.slopes)))]  # replace remaining empty values (from the first few days of observations) with the value from the first date where I have adjustment values

plot(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0] ~ as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0]))))  # plot the calculated date-specific slopes
lines(min(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])))):max(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])))), predict(slope.by.date.loess, min(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0])))):max(as.numeric(as.Date(names(caliper.adjust$date.slopes[caliper.adjust$date.slopes!=0]))))))  # plot the loess curve

rm(list=c("caliper.variation.optimized_slope.averaged", grep("^dates", x=ls(), value=TRUE), grep("^slope", x=ls(), value=TRUE), grep("^inst_catnum", x=ls(), value=TRUE))) # removed unneeded variables from this section


## finally, adjust the data by the date-specific adjustments
picidae.morph.log.reduced_var$adjusted_date <- picidae.morph.log.reduced_var$base  # create a copy of the data, to perform adjustments on this copy
for (i in 1:nrow(picidae.morph.log.reduced_var$adjusted_date)) {  # loop over all sets of measurement data (rows in the data frame)
  picidae.morph.log.reduced_var$adjusted_date[i,-(1:5)] <- (picidae.morph.log.reduced_var$adjusted_date[i,-(1:5)] - caliper.adjust$intercept) / caliper.adjust$date.slopes[as.character(picidae.morph.log.reduced_var$adjusted_date$Data_date[i])]
}  # adjust the data for the current row, using the intercept and the date-specific slope
rm(i)


### take means by individual (includes forking data by individual completeness, complete_ind_only vs. all), merging caliper data with data from images, and taking means by sex and taxon (includes forking data by taxon completeness, mf_only vs. all - mf_only includes only taxa for which I have both males AND females)

## calculate the mean of all measurements for each individual (for cases where the same individual had the same character measured >1 times)
picidae.morph.log.reduced_var_inds <- list()  # create a list to store data reduced to a single row per individual
picidae.morph.log.reduced_var_inds$adjusted_date <- reduceToMeanByCols(picidae.morph.log.reduced_var$adjusted_date, reduce_colnames=c("Inst_CatNum","My_genus_species","Sex"), mean_cols=-(1:5))  # take means by individual
picidae.morph_from_images$log.reduced_inds$Inst_CatNum[!(picidae.morph_from_images$log.reduced_inds$Inst_CatNum %in% picidae.morph.log.reduced_var_inds$adjusted_date$Inst_CatNum)]  # check for individuals in the measurements from images but not in the measurements from calipers; this avoids picking up those that were excluded earlier; there should be few to none of these, and might indicate mismatched Inst_CatNum between the image data and caliper data
picidae.morph.log.reduced_var_inds$adjusted_date$Inst_CatNum[!(picidae.morph.log.reduced_var_inds$adjusted_date$Inst_CatNum %in% picidae.morph_from_images$log.reduced_inds$Inst_CatNum)]  # check for individuals in the measurements from calipers but not in the measurements from images; this avoids picking up those that were excluded earlier; these should be individuals that had no skulls or unmeasurable skulls

## combine the data from caliper measurements (reduced for variables and individuals, and adjusted for date variation) with the data from image measurements, keeping all of the caliper ones (even if they don't have corresponding image measurements)
picidae.morph_combined.log.reduced_var_inds <- list()  # create a list to store combined image and caliper data
picidae.morph_combined.log.reduced_var_inds$all_inds <- merge(picidae.morph.log.reduced_var_inds$adjusted_date, subset(picidae.morph_from_images$log.reduced_inds, Inst_CatNum %in% picidae.morph.log.reduced_var_inds$adjusted_date$Inst_CatNum), by="Inst_CatNum", all.x = TRUE)  # combine the data from caliper measurements (reduced for variables and individuals, and adjusted for date variation) with the data from image measurements, keeping all of the caliper ones (even if they don't have corresponding image measurements)
picidae.morph_combined.log.reduced_var_inds$all_inds <- cbind(picidae.morph_combined.log.reduced_var_inds$all_inds[1:3], picidae.morph_combined.log.reduced_var_inds$all_inds[grep("ramp", names(picidae.morph_combined.log.reduced_var_inds$all_inds))], picidae.morph_combined.log.reduced_var_inds$all_inds[-c(1:3, grep("ramp", names(picidae.morph_combined.log.reduced_var_inds$all_inds)))]) # re-order the columns in the data frame

# generate a separate data frame containing only individuals with complete data; considered complete if they have either ramp_absent or ramp_present for the various skull measurements and either left or right for the post-cranial measurements
picidae.morph_combined.log.reduced_var_inds$complete_ind_only <- subset(picidae.morph_combined.log.reduced_var_inds$all_inds, subset=((!is.na(Cranium_GL_ramp_present) | !is.na(Cranium_GL_ramp_absent)) & (!is.na(Cranium_OC_ramp_present) | !is.na(Cranium_OC_ramp_absent)) & (!is.na(Mandible_mean_ramp_present) | !is.na(Mandible_mean_ramp_absent)) & !is.na(Humerus_mean) & !is.na(Ulna_mean) & !is.na(Femur_mean) & !is.na(Tibiotarsus_mean) & !is.na(Tarsometatarsus_mean) & !is.na(Pygostyle_width) & !is.na(Pygostyle_length) & (!is.na(Maxilla_length_ramp_absent) | !is.na(Maxilla_length_ramp_present))))

## find the mean of all measurements for each sex for each species; basically, this reduces the data to one average representative per sex per species
picidae.morph_combined.log.reduced_var_inds_sex <- list() # create a list to store the data reduced to one row per sex per species
for (i in names(picidae.morph_combined.log.reduced_var_inds)) {  # loop over individual inclusion
    picidae.morph_combined.log.reduced_var_inds_sex[[i]] <- reduceToMeanByCols(picidae.morph_combined.log.reduced_var_inds[[i]], reduce_colnames=c("My_genus_species","Sex"), mean_cols=4:ncol(picidae.morph_combined.log.reduced_var_inds[[i]]))
}
rm(i)

## summarize data frame to determine if a species has both M and F; if so, just average those; if not, include all individuals in the average
table.picidae.morph_combined.taxa_by_sex <- list()  # create a list to store tables of sex by species
for (i in names(picidae.morph_combined.log.reduced_var_inds_sex)) {  # loop over individual inclusion
    table.picidae.morph_combined.taxa_by_sex[[i]] <- table(picidae.morph_combined.log.reduced_var_inds_sex[[i]][,1:2])  # generate table of taxon (row) by sex (column), so that I can check which taxa have males and females
}
rm(i)

## remove taxa that I won't be using in the analyses, as I don't want them to contribute to the data interpolation
picidae.taxa.exclude <- list()  # create a list to store taxa to exclude
picidae.taxa.exclude$not_species <- c("Picidae_misidentified", grep(pattern="/|_x_", rownames(table.picidae.morph_combined.taxa_by_sex$all_vars$all_inds), value=TRUE)) # generate a vector of ambiguous taxa to remove from the morphological data set; this finds hybrids and ambiguous taxa based on the "/" and "_x_" patterns in the names
picidae.taxa.exclude$no_phylo <- setdiff(unique(picidae.morph$log$My_genus_species)[order(unique(picidae.morph$log$My_genus_species))], c(picidae.RAxML.all.BEAST_calibrated$tip.label, picidae.taxa.exclude$not_species))  # identify taxa with data but not in phylogeny (and not non-species taxa); treedata does this, but I need it in general for possible use of taxonomic proxies
picidae.taxa.exclude$no_data <- setdiff(picidae.RAxML.all.BEAST_calibrated$tip.label, unique(picidae.morph$log$My_genus_species))  # identify taxa in the phylogeny but with no measurement data


### checking if there are species I should use as proxies, where I have phylogenetic data for only one member of a clear species pair or group, and measurement data for only one (but different) member of that same species pair or group
# if these are present, treat them as if the phylogenetic information were for the taxon for which I have morph data; change the name in the tree (but make a note, and make a copy of the tree that's different)

# replace tip names for a few taxa in the phylogeny so that I can include one member of a species pair or species group (where I have phylo for one and data for another)
picidae.RAxML.all.BEAST_calibrated.with_proxies <- picidae.RAxML.all.BEAST_calibrated
picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label[picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label=="Gecinulus_grantia"] <- "Gecinulus_viridis"
picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label[picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label=="Dendropicos_pyrrhogaster"] <- "Dendropicos_namaquus"
picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label[picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label=="Melanerpes_herminieri"] <- "Melanerpes_portoricensis"
picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label[picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label=="Piculus_callopterus"] <- "Piculus_simplex"

## extract lists of taxa with data for males and females, and taxa missing data for either males or females
picidae.morph_combined.taxa_by_mf <- list() # generate a list to store taxon lists for data variants
for (i in names(table.picidae.morph_combined.taxa_by_sex)) {  # loop over individual inclusion
    picidae.morph_combined.taxa_by_mf[[i]][["no_female_male"]] <- rownames(table.picidae.morph_combined.taxa_by_sex[[i]][table.picidae.morph_combined.taxa_by_sex[[i]][,"F"]==0 | table.picidae.morph_combined.taxa_by_sex[[i]][,"M"]==0,]) # generate a vector of taxa that are missing either males OR females OR both
    picidae.morph_combined.taxa_by_mf[[i]][["female_male"]] <- rownames(table.picidae.morph_combined.taxa_by_sex[[i]][table.picidae.morph_combined.taxa_by_sex[[i]][,"F"]==1 & table.picidae.morph_combined.taxa_by_sex[[i]][,"M"]==1,]) # generate a vector of taxa that have both males AND females
}
rm(i)


## find the mean of all measurements for each taxon, averaging across sexes
picidae.morph_combined.log.reduced_var_inds_sex_taxa <- list()
for (i in names(picidae.morph_combined.taxa_by_mf)) {  # loop over individual inclusion
  picidae.morph_combined.log.reduced_var_inds_sex_taxa[[i]] <- reduceToMeanByCols(data = subset(picidae.morph_combined.log.reduced_var_inds_sex[[i]], subset=((Sex %in% c("F","M")) & (My_genus_species %in% picidae.morph_combined.taxa_by_mf[[i]]$female_male) & !(My_genus_species %in% picidae.taxa.exclude$not_species))), reduce_colnames="My_genus_species", mean_cols=3:ncol(picidae.morph_combined.log.reduced_var_inds_sex[[i]]))  # calculate means by species with only species where I have both males and females; it checks against the vector containing taxa with both males and females; it also only includes rows of data where the sex is "M" or "F"
  picidae.morph_combined.log.reduced_var_inds_sex_taxa[[i]] <- rbind(picidae.morph_combined.log.reduced_var_inds_sex_taxa[[i]], reduceToMeanByCols(data = subset(picidae.morph_combined.log.reduced_var_inds[[i]], subset=(My_genus_species %in% picidae.morph_combined.taxa_by_mf[[i]]$no_female_male) & !(My_genus_species %in% picidae.taxa.exclude$not_species)), reduce_colnames="My_genus_species", mean_cols=4:ncol(picidae.morph_combined.log.reduced_var_inds[[i]])))  # add in the species that are missing either males or females or both, using the averages of ALL INDIVIDUALS, not just the mean of each sex
}
rm(i)


### impute data without ramphotheca from measurements with ramphotheca (to standardize the measurements)

## the function imputeRampAbsent() fits models of ramp_absent data to ramp_present data, and imputes missing data (for species or individuals) using stochastic imputation based on those models (with a number of options)
# as input, it takes data_by_species (the data with means by species), impute_fields (names of fields to impute), data_by_inds (the data with means by individual, unnecessary if species_inds="species"), model_fitting ("pgls" for PGLS model, "lm" for non-phylogenetic model), phy (the phylogeny, required if fitting with PGLS), species_inds ("species" to impute by species, "inds" to impute by individuals), return_models (boolean for whether to return the results of model fitting), return_treedata (boolean for whether to return the treedata objects used in PGLS fitting)
# it returns a list containing a data frame with imputed data where it was absent for ramp_absent, imputed either by individuals or by species, optionally the model fits and values of pseudo-R^2, and optionally treedata objects used in PGLS mdoel fitting
imputeRampAbsent <- function(data_by_species, impute_fields, data_by_inds=NULL, model_fitting="pgls", phy=NULL, species_inds="species", return_models=FALSE, return_treedata=FALSE) {
  require(ape)
  require(nlme)
  require(geiger)
  require(stringr)
  
  # make a copy of the data frame to add imputed values to
  if (species_inds=="species") {data.imputed <- data_by_species
  } else if (species_inds=="inds") {data.imputed <- data_by_inds}
  
  models = list()  # initialize a list to store model results
  treedata_by_field = list()  # initialize a list to store treedata objects
  cors = numeric()  # initialize a vector to store pesudo-R^2 values
  data_by_species.matrix <- data.matrix(data_by_species[,-1]) # create a matrix version of the original species-level data (in order to use geiger::treedata)
  rownames(data_by_species.matrix) <- data_by_species$My_genus_species # set the row names of the matrix to be the taxon names (treedata needs this)
  
  for (i in 1:length(impute_fields)) {  # loop over fields to impute
    field.ramp_absent <- paste(impute_fields[i], "_ramp_absent", sep="")
    field.ramp_present <- paste(impute_fields[i], "_ramp_present", sep="")
    indices.to_fill.tmp <- which(is.na(data.imputed[,field.ramp_absent]) & !is.na(data.imputed[,field.ramp_present])) # extract the indices of the data that need to be imputed (so that they don't need to be re-evaluated them every time, and so that they can be used even after the values are filled int (e.g. to add residuals))
    
    if (model_fitting=="pgls") {
      data.tmp.treedata <- treedata(phy=phy, data=data_by_species.matrix[!is.na(data_by_species.matrix[,field.ramp_absent]) & !is.na(data_by_species.matrix[,field.ramp_present]), c(field.ramp_absent, field.ramp_present)], sort=TRUE, warnings=FALSE) # create a treedata object that only has taxa and data for which the current field has data for ramp_present and ramp_absent
      model.tmp <- gls(data = data.frame(y = data.tmp.treedata$data[,field.ramp_absent], x = data.tmp.treedata$data[,field.ramp_present]), model = y ~ x, na.action=na.exclude, correlation=corPagel(value=1, phy=data.tmp.treedata$phy), method="REML")  # fit PGLS model, with lambda estimated
      
    } else if (model_fitting=="lm") {
      model.tmp <- lm(data = data.frame(y = data_by_species.matrix[!is.na(data_by_species.matrix[,field.ramp_absent]) & !is.na(data_by_species.matrix[,field.ramp_present]),field.ramp_absent], x = data_by_species.matrix[!is.na(data_by_species.matrix[,field.ramp_absent]) & !is.na(data_by_species.matrix[,field.ramp_present]),field.ramp_present]), formula = y ~ x, na.action=na.exclude) # fit simple linear model assuming independent data
    }
    
    data.imputed[indices.to_fill.tmp, field.ramp_absent] <- predict(model.tmp, newdata=data.frame(x = data.imputed[indices.to_fill.tmp, field.ramp_present]), na.action=na.pass) # find the values that are NA for the current field with ramp_absent, and replace them (where possible) with predicted values from the current model (still need to add stochastic variation)
    
    # calculate the pseudo-residuals for species not included in pgls models due to lack of phylogenetic information, and append them to the model residuals (this allows me to have them incorporated downstream without any changes)
    if (model_fitting=="pgls") {
      species.resids.no_phylo <- data_by_species$My_genus_species[!(data_by_species$My_genus_species %in% names(model.tmp$resid)) & (!is.na(data_by_species[[field.ramp_absent]])) & (!is.na(data_by_species[[field.ramp_present]]))]  # identify species with data for both ramp_absent and ramp_present but not in the phylogeny
      resids.no_phylo <- data_by_species[[field.ramp_absent]][match(species.resids.no_phylo, data_by_species$My_genus_species)] - (model.tmp$coefficients[2] * data_by_species[[field.ramp_present]][match(species.resids.no_phylo, data_by_species$My_genus_species)] + model.tmp$coefficients[1])  # determine species-level pseudo-residuals using empirical values and predicted values from model
      names(resids.no_phylo) <- species.resids.no_phylo  # assign species names to calculated pseudo-residuals
      model.tmp$resid <- c(model.tmp$resid, resids.no_phylo)  # add calculated pseudo-residuals to vector of model residuals
    }
    
    # add residuals to individual-level or species-level data using species-level residuals, if present, then using the genus average residuals, if the species-level one isn't present
    if (species_inds=="inds") {
      indices.to_fill.sp.resids <- indices.to_fill.tmp[!is.na(model.tmp$resid[data.imputed$My_genus_species[indices.to_fill.tmp]])]  # get the indices for which I have species-level residuals
      data.imputed[indices.to_fill.sp.resids, field.ramp_absent] <- data.imputed[indices.to_fill.sp.resids, field.ramp_absent] + model.tmp$resid[data.imputed$My_genus_species[indices.to_fill.sp.resids]] + rnorm(sd=summary(model.tmp)$sigma)  # add the species-level residuals to the data, and a random variable using the residual standard error from the model
      indices.to_fill.tmp <- indices.to_fill.tmp[is.na(model.tmp$resid[data.imputed$My_genus_species[indices.to_fill.tmp]])]  # modify the vector of indices so that it only contains those without species-level data
    }
    genera.to_fill.tmp <- str_extract(data.imputed[indices.to_fill.tmp, "My_genus_species"], pattern="[A-Za-z]+")  # pull out the genus names for the data to fill
    resid.genera <- str_extract(names(model.tmp$resid), pattern="[A-Za-z]+")  # get the genus names from the model residuals (to facilitate matching)
    indices.to_fill.genus.resids <- indices.to_fill.tmp[genera.to_fill.tmp %in% resid.genera]  # cut the vector of indices down to only include those for which I have residuals from other genus members
    genera.to_fill.resids <- genera.to_fill.tmp[genera.to_fill.tmp %in% resid.genera]  # get the genera that need to have residuals added
    for (j in 1:length(indices.to_fill.genus.resids)) {  # cut the vector of genera down to only include those for which I have residuals from other genus members
      data.imputed[indices.to_fill.genus.resids[j], field.ramp_absent] <- data.imputed[indices.to_fill.genus.resids[j], field.ramp_absent] + mean(model.tmp$resid[resid.genera == genera.to_fill.resids[j]]) + rnorm(sd=summary(model.tmp)$sigma)  # add the genus average residuals (where species-level residuals not available), and a random variable using the residual standard error from the model
    }
    
    # store additional data to return (optional)
    if (return_models) {
      models[[i]] <- model.tmp # store the model
      if (model_fitting=="pgls") {
        cors[i] <- cor(data.tmp.treedata$data[,1],predict(model.tmp))^2 # calculate a pseudo-R2 for this fit
      } else if (model_fitting=="lm") {
        cors[i] <- summary(model.tmp)$r.squared
      }
      names(models)[i] <- names(cors)[i] <- impute_fields[i]
    }
    if (return_treedata & model_fitting=="pgls") {
      treedata_by_field <- data.tmp.treedata
      names(treedata_by_field)[i] <- impute_fields[i]
    }
  }
  
  data.imputed <- data.imputed[,-which(colnames(data.imputed)==field.ramp_present)] # drop ramp_present columns from data frame
  
  # return the appropriate objects
  if (!(return_treedata & model_fitting=="pgls")) {
    if (!return_models) {
      return(list(data.imputed = data.imputed))
    } else {
      return(list(data.imputed = data.imputed, models = models, cors = cors))
    }
  } else {
    if (!return_models) {
      return(list(data.imputed = data.imputed, treedata=treedata_by_field))
    } else {
      return(list(data.imputed = data.imputed, models = models, cors = cors, treedata=treedata_by_field))
    }
  }
}

## generate new data frames to store imputed values for individuals with ramphotheca absent
picidae.morph_combined.log.reduced_var_inds.imputed <- list()
for (i in names(picidae.morph_combined.log.reduced_var_inds_sex_taxa)) {  # loop over individual inclusion
    picidae.morph_combined.log.reduced_var_inds.imputed[[i]] <- imputeRampAbsent(data_by_species=picidae.morph_combined.log.reduced_var_inds_sex_taxa[[i]], data_by_inds=picidae.morph_combined.log.reduced_var_inds[[i]], phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, impute_fields = c("Cranium_GL", "Cranium_OC", "Mandible_mean", "Maxilla_length"), model_fitting="pgls", species_inds="inds", return_models=TRUE, use_resids=l)
}
rm(i)

## now re-aggregate the data by sex, exactly as I did it before
picidae.morph_combined.log.reduced_var_inds_sex.imputed <- list()
for (i in names(picidae.morph_combined.log.reduced_var_inds.imputed)) {  # loop over individual inclusion
    picidae.morph_combined.log.reduced_var_inds_sex.imputed[[i]] <- reduceToMeanByCols(picidae.morph_combined.log.reduced_var_inds.imputed[[i]][["data.imputed"]], reduce_colnames=c("My_genus_species", "Sex"), mean_cols=4:ncol(picidae.morph_combined.log.reduced_var_inds.imputed[[i]][["data.imputed"]]))
}
rm(i)

## now re-aggregate the data by taxon, exactly as I did it before
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed <- list()
for (i in names(picidae.morph_combined.log.reduced_var_inds_sex.imputed)) {  # loop over individual inclusion
  picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]] <- reduceToMeanByCols(data = subset(picidae.morph_combined.log.reduced_var_inds_sex.imputed[[i]], subset=((Sex %in% c("F","M")) & (My_genus_species %in% picidae.morph_combined.taxa_by_mf[[i]][["female_male"]]) & !(My_genus_species %in% picidae.taxa.exclude$not_species))), reduce_colnames="My_genus_species", mean_cols = 3:ncol(picidae.morph_combined.log.reduced_var_inds_sex.imputed[[i]]))  # calculate means by species with only species where I have both males and females; it checks against the vector containing taxa with both males and females; it also only includes rows of data where the sex is "M" or "F"
  picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]] <- rbind(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]], reduceToMeanByCols(data = subset(picidae.morph_combined.log.reduced_var_inds.imputed[[i]][["data.imputed"]], subset=(My_genus_species %in% picidae.morph_combined.taxa_by_mf[[i]][["no_female_male"]]) & !(My_genus_species %in% picidae.taxa.exclude$not_species)), reduce_colnames="My_genus_species", mean_cols=4:ncol(picidae.morph_combined.log.reduced_var_inds.imputed[[i]][["data.imputed"]])))
  picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]] <- picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]][order(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]]$My_genus_species),] # sort by My_genus_species, to make things simpler
}
rm(i)


### check data again (after interpolation, in part to check interpolation)

k <- "all_inds"

## plot histograms of all measurements for each species, to check data for outliers
for (i in picidae.species.names) {
  pdf(file=paste("/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Picidae_combined_analyses_dist_morph_evol_divers/Picidae_morph_data_checking/interpolated_data_checking/", gsub("/", "_", i), "_histogram_matrix.pdf", sep=""), width=10, height=8)
  par(mfrow=c(3,4))
  for (j in 4:14) {
    data.tmp <- picidae.morph_combined.log.reduced_var_inds.imputed[[k]][["data.imputed"]][picidae.morph_combined.log.reduced_var_inds.imputed[[k]][["data.imputed"]]$My_genus_species==i, j]
    if (any(!is.na(data.tmp))) {hist(data.tmp, main=colnames(picidae.morph_combined.log.reduced_var_inds.imputed[[k]][["data.imputed"]])[j], breaks=10)
    } else plot.new()
  }
  dev.off()
}
rm(i,j,k,data.tmp)


### convert to complete data set (remove taxa with not enough data, and impute data for those with only a single measurement missing)

findTaxaExclude <- function(data, threshold=1, ident="My_genus_species") {
  # this function counts up the NAs in a data frame, and returns a vector of the identifiers (as ident) of members that have more NAs than threshold
  na_counts <- rep(NA, (nrow(data)))
  for (i in 1:nrow(data)) {
    na_counts[i] <- sum(is.na(data[i,]))
  }
  print(paste(data[na_counts > 0, ident], " : ", na_counts[na_counts > 0], sep=""))
  taxa_to_exclude <- data[na_counts > threshold, ident]
  return(taxa_to_exclude)
}

# identify taxa that don't have enough data to use
picidae.taxa.exclude$not_enough_data  <- findTaxaExclude(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[["all_inds"]])  # taxa to remove due to not enough data (anything with NA in more than one field); this is the same regardless of how I impute data or treat residuals, and there are none for complete_ind

## create a new list of data frames, minus the taxa that don't have enough data
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa <- list()
for (i in names(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed)) {  # loop over individual inclusion
  picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa[[i]] <-  picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa[[i]] <- subset(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]], !(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[i]]$My_genus_species %in% picidae.taxa.exclude$not_enough_data[[i]])) # modify the data frame to remove those taxa; no need to do this for complete_ind_only, as all taxa have all data in that set, but I did it for all, so that I would have one consistent data object to use
}
rm(i)


### impute missing data, where one measurement is missing from a given taxon

## the function imputeMissingVars() uses stochastic imputation to impute data for species missing data for only one variable; it loops over each variable that has NAs, building models from all the other variables (complete cases only), then filling in the NAs based on those models where possible
# as input, it takes data (as a data frame containing species means), phy (phylogeny, required if model_fitting="pgls"), model_fitting (if "pgls" fits a PGLS model, if "lm" fits a non-phylogenetic model), return_models (boolean for whether to return model fits)
imputeMissingVars <- function(data, phy=NULL, model_fitting="pgls", return_models=FALSE) {
  require(ape)
  require(geiger)
  require(nlme)
  
  data.imputed <- data  # create a copy of the data for imputing
  
  # if there are no NAs in the data, skip the rest and return the original data frame
  if (sum(is.na(data)) == 0) {
    paste("No NAs found; returning original data frame")
    return(list(data.imputed = data.imputed))
  }
  
  impute_fields <- colnames(data)[which(colSums(is.na(data)) > 0)] # find the data fields to interpolate
  models = list()  # initialize a list to store the models
  cors = numeric()  # initialize a vector to store the pseudo-R^2 values
  
  if (model_fitting=="pgls") {
    data.matrix.na_omitted <- data.matrix(data[which(rowSums(is.na(data)) == 0),-1]) # create a matrix version of the data with only complete cases (in order to use geiger::treedata); omit NAs because I'm only building the models from complete cases, and any NA means an incomplete case
    rownames(data.matrix.na_omitted) <- na.omit(data)$My_genus_species  # set rownames as the species (necessary for treedata)
    data.treedata <- treedata(phy=phy, data=data.matrix.na_omitted, sort=TRUE, warnings=FALSE)  # create a treedata object
    
  } else if (model_fitting=="lm") {
    data.na_omitted <- data[which(rowSums(is.na(data)) == 0),-1]  # create a matrix version of the data with only complete cases
  }
  
  for (i in 1:length(impute_fields)) {  # loop over fields to impute
    indices.to_fill.tmp <- which(is.na(data[,impute_fields[i]]))  # identify indices of missing data to impute
    
    if (model_fitting=="pgls") {
      data.tmp.data_frame <- data.frame(data.treedata$data)[,colnames(data.treedata$data) != impute_fields[i]]  # create a data frame to store data for the current model (makes model fitting easier)
      data.tmp.data_frame$y <- data.treedata$data[,impute_fields[i]]  # add respone variable to data frame
      model.tmp <- gls(data = data.tmp.data_frame, model = y~., na.action=na.exclude, correlation=corPagel(value=1, phy=data.treedata$phy), method="REML") # fit PGLS model to the data, regressing the response variable on all other variables, with lambda estimated
    } else if (model_fitting=="lm") {
      data.tmp.data_frame <- data.na_omitted[, colnames(data.na_omitted) != impute_fields[i]]    # create a data frame to store data for the current model (makes model fiting easier)
      data.tmp.data_frame$y <- data.na_omitted[, impute_fields[i]]  # add respone variable to data frame
      model.tmp <- lm(data = data.tmp.data_frame, formula = y~.)  # fit linear model to the data, regressing the response variable on all other variables
    }
    
    data.imputed[indices.to_fill.tmp, impute_fields[i]] <- predict(model.tmp, newdata=data[indices.to_fill.tmp, !(colnames(data.imputed) %in% c(impute_fields[i], "My_genus_species"))]) + rnorm(n=length(indices.to_fill.tmp), sd=summary(model.tmp)$sigma)  # impute data, using the predicted value based on all other variables, and a random variable using the residual standard error from the model 
    
    if (return_models) {
      models[[i]] <- model.tmp # store the model
      if (model_fitting=="pgls") {
        cors[i] <- cor(data.treedata$data[,impute_fields[i]],predict(model.tmp))^2 # calculate a pseudo-R2 for this fit
      } else if (model_fitting=="lm") {
        cors[i] <- summary(model.tmp)$r.squared
      }
      names(models)[i] <- names(cors)[i] <- impute_fields[i]
    }
  }
  
  if (!return_models) {
    return(list(data.imputed = data.imputed))
  } else {
    return(list(data.imputed = data.imputed, models = models, cors = cors))
  }
}

## create new data frame with data interpolated where only one measurement is missing
# don't need to do this for complete_ind_only, because they don't have any NAs; but I'm including those in this list so that I have one consistent data object to use
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa.imputed <- list()
for (i in names(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa)) {  # loop over individual inclusion
    picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa.imputed[[i]] <- imputeMissingVars(data=picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa[[i]], phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, model_fitting="pgls", return_models=TRUE)
}
rm(i)

picidae.morph.log.fully_reduced <- picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa.imputed  # store a copy of the fully reduced data set with a simplified name


### create data variants (geomean, treedata objects, size-scaled shape variables)

## the function sizeScalePGLS() takes a data frame of log-transformed data values, with values to scale in scale_cols, and a tree containing the taxa those data represent, and adds a geomean column to the end, then scales the variables in scale_cols by the geomean (by finding the residuals from PGLS regression of each variable on the geomean)
#  it can do this with incomplete overlap in taxa; it uses treedata to drop the taxa not found in both, and outputs a list of the taxa that were dropped
# as input, it takes data (a data frame containing species means), scale_cols (vector of numbers of columns to scale), phy (phylogeny), name_col (name or number of the column of the data frame with the taxon names), keep_geomean (boolean to retain the geomean column in the returned data frame), return_models (boolean to return the model fits), and optional methods to pass to the gls function for PGLS model fitting
# it returns either a data frame with the size-scaled data (with or without geomean) or a list containing the size-scaled data and the model fits
sizeScalePGLS <- function(data, scale_cols, phy, name_col, keep_geomean=FALSE, return_models=FALSE, method="ML", control=list()) {
  require(geiger)
  require(nlme)
  scale_cols <- scale_cols[order(scale_cols)]  # put the columns to scale in order (makes things easier)
  data$geomean <- apply(data[,scale_cols],MARGIN=1,mean)  # add a column to the data frame, which is the arithmetic mean of the log-transformed values (which is the same as the log of the geometric mean of the untransformed values)
  data.as_matrix <- data.matrix(data[,c(scale_cols, which(colnames(data)=="geomean"))])  # create a matrix version of the data (in order to use geiger::treedata)
  rownames(data.as_matrix) <- data[,name_col]  # add the taxon names as rownames
  data.treedata <- treedata(phy=phy, data=data.as_matrix, sort=TRUE)  # create a treedata object to filter out taxa missing from the data or the phylogeny
  data.geomean_scaled <- data[match(data.treedata$phy$tip.label, data[,name_col]),]  # create a copy of the data frame, including only those taxa which are included in the treedata object (and therefore are present in both the phy and the data), in the order of the tip.labels, so that the order is consistent
  
  # fit PGLS models for each variable regressed on geomean; it uses try() in order to catch occasional problems with fitting the lambda model
  model.pgls.traits_geomean <- list()  # initialize a list to store model fits
  for (i in 1:length(scale_cols)) {  # loop over columns to scale
    model.tmp <- try(gls(data=data.frame(x = data.treedata$data[,i], geomean = data.treedata$data[,"geomean"]), model = x ~ geomean, correlation=corPagel(value=1, phy=data.treedata$phy), method=method, control=control), silent=TRUE)  # fit a PGLS model with lambda estimated
    if (class(model.tmp) == "gls") {
      model.pgls.traits_geomean[[i]] <- model.tmp  # if the fit worked, pass it on
    } else {
      cat("corPagel failed on column ", colnames(data)[scale_cols[i]], ". Trying corBrownian.\n", sep="")
      cat("Error: ", model.tmp, "\n", sep="")
      model.tmp <- try(gls(data=data.frame(x = data.treedata$data[,i], geomean = data.treedata$data[,"geomean"]), model = x ~ geomean, correlation=corBrownian(value=1, phy=data.treedata$phy), method=method, control=control), silent=TRUE)  # if the fit didn't work, try using a Brownian motion model (instead of estimating lambda)
      if (class(model.tmp) == "gls") {
        model.pgls.traits_geomean[[i]] <- model.tmp  # if the fit worked, pass it on
      } else {
        # if the fit didn't work, stop
        cat("corPagel and corBrownian both failed. Try something else.\n")
        cat("Error: ", model.tmp, "\n", sep="")
        return()
      }
    }
  }
  
  # replace the original values with the residuals from the PGLS regression of each variable on geomean
  for (i in 1:length(scale_cols)) {  # loop over columns to scale
    data.geomean_scaled[,scale_cols[i]] <- model.pgls.traits_geomean[[i]]$resid  # replace values with residuals from the model fitting
  }
  rm(i)
  
  if (!keep_geomean) data.geomean_scaled$geomean <- NULL  # remove geomean column if specified
  
  if (return_models) {
    names(model.pgls.traits_geomean) <- names(data)[scale_cols]
    return(list(data.geomean_scaled=data.geomean_scaled, model.pgls.traits_geomean=model.pgls.traits_geomean))
  } else return(data.geomean_scaled)
}


## create treedata objects for each of the data sets (mostly so that I have an easy-to-access tree for each data variant)
picidae.morph.log.fully_reduced.treedata <- list()
for (i in names(picidae.morph.log.fully_reduced)) {  # loop over individual inclusion
  data.matrix.tmp <- data.matrix(picidae.morph.log.fully_reduced[[i]][["data.imputed"]][,2:ncol(picidae.morph.log.fully_reduced[[i]][["data.imputed"]])])
  rownames(data.matrix.tmp) <- picidae.morph.log.fully_reduced[[i]][["data.imputed"]]$My_genus_species
  picidae.morph.log.fully_reduced.treedata[[i]] <- treedata(data=data.matrix.tmp, phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, sort=TRUE)
}
rm(i,data.matrix.tmp)

## create a separate object with the geomean
picidae.morph.log.fully_reduced.geomean <- list()
for (i in names(picidae.morph.log.fully_reduced)) {  # loop over individual inclusion
  picidae.morph.log.fully_reduced.geomean[[i]] <- apply(picidae.morph.log.fully_reduced[[i]][["data.imputed"]][,2:ncol(picidae.morph.log.fully_reduced[[i]][["data.imputed"]])], MARGIN=1, mean)
  names(picidae.morph.log.fully_reduced.geomean[[i]]) <- picidae.morph.log.fully_reduced[[i]][["data.imputed"]]$My_genus_species
}
rm(i)

## create an additional morphology data frame, with a column for the geometric mean (across all measurements within each taxon), and all other measurements scaled by the geometric mean of that taxon (using residuals from PGLS regression of each variable on the geomean)
picidae.morph.log.fully_reduced.geomean_scaled <- list()
for (i in names(picidae.morph.log.fully_reduced)) {  # loop over individual inclusion
  cat("Starting", i, "\n", sep=" ")
  picidae.morph.log.fully_reduced.geomean_scaled[[i]] <- sizeScalePGLS(data = picidae.morph.log.fully_reduced[[i]][["data.imputed"]], scale_cols=2:ncol(picidae.morph.log.fully_reduced[[i]][["data.imputed"]]), phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, name_col=1, return_models=TRUE, method="ML")
}
rm(i,j,k,l)

## plot them all against each other to see correlations; these scatterplot matrices are super interesting
for (i in names(picidae.morph.log.fully_reduced.geomean_scaled)) {  # loop over individual inclusion
  plot(picidae.morph.log.fully_reduced.geomean_scaled[[i]][["data.geomean_scaled"]][,-1], main=i) # plot them all against each other to see correlations
}
rm(i)


### run phylogenetic PCA on the morphological data

## the function wrappedPhylPCA() is a wrapper to run phylogenetic PCA; it creates the necessary intermediate steps (matrix, treedata), and also calculates the percentage of the variance captured by each PCA
# as input, it takes data (the data as a data frame), pca_cols (the columns on which to conduct phylogenetic PCA), phy (the phylogeny), name_col (the name or number of the column containing taxon names), and method (the method to pass to phyl.pca)
# it returns a list containing the results of the phylogenetic PCA (an object of class phyl.pca) and a vector of the percent variance explained by each PCA axis
wrappedPhylPCA <- function(data, pca_cols, phy, name_col, method="lambda") {
  require(geiger)
  require(phytools)
  newdata.matrix <- data.matrix(data[,pca_cols])  # create a matrix version of the data (in order to use geiger::treedata)
  rownames(newdata.matrix) <- data[,name_col]  # set the rownames to be the taxon names (treedata needs this)
  newdata.treedata <- treedata(phy = phy, data = newdata.matrix, warnings=!quiet, sort=TRUE)  # construct a treedata object (which reduces the phylogeny and the data matrix based on shared membership)
  newdata.pca <- phyl.pca(tree = newdata.treedata$phy, Y = newdata.treedata$data, method=method, mode="cov")  # generate a phylogenetic PCA of the values
  newdata.variances <- diag(newdata.pca$Eval)/sum(newdata.pca$Eval)*100  # calculate the percentage of the variance in the data that is captured be each PCA axis
  return(list(pca = newdata.pca, var_percents = newdata.variances))
}


## run phylogenetic PCA on the unscaled morphological data
picidae.morph.log.fully_reduced.phyl_pca <- list()
for (i in names(picidae.morph.log.fully_reduced)) {  # loop over individual inclusion
  cat("Starting", i, "\n", sep=" ")
  picidae.morph.log.fully_reduced.phyl_pca[[i]] <- wrappedPhylPCA(data = picidae.morph.log.fully_reduced[[i]][["data.imputed"]], pca_cols=2:ncol(picidae.morph.log.fully_reduced[[i]][["data.imputed"]]), phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, name_col=1)  # conduct the phylogenetic on the current data set
}
rm(i)

## run phylogenetic PCA on the size-scaled shape morphological data
picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca <- list()
for (i in names(picidae.morph.log.fully_reduced.geomean_scaled)) {  # loop over individual inclusion
  cat("Starting", i, "\n", sep=" ")
  picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]] <- wrappedPhylPCA(data = picidae.morph.log.fully_reduced.geomean_scaled[[i]][["data.geomean_scaled"]], pca_cols=2:ncol(picidae.morph.log.fully_reduced.geomean_scaled[[i]][["data.geomean_scaled"]]), phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, name_col=1)  # conduct the phylogenetic on the current data set
}
rm(i)


### generate and output trees for use in downstream analyses

## output the tree containing only the taxa for which I have phylogenetic and morphological data (for use in other analyses e.g. BAMM)
write.tree(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy, file="RAxML_bestTree_BEAST_calibrated_morphdata_only.tre")

## plot the trees containing only the taxa for which I have phylogenetic and morphological data
pdf(file="RAxML_bestTree_BEAST_calibrated_morphdata_only.pdf", width=4.5, height=10)
plot(ladderize(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy, right=TRUE), direction="rightwards", cex=0.45, no.margin=T, x.lim=c(-2, 45), label.offset=0.2)
dev.off()


### PCA-rotate and geomean-scale the data without phylo, using the rotations and/or models built from data with phylo


## the function sizeScalePGLSNoPhylo() uses the models built from PGLS regression of measurements on size to scale data not included in those models because the species are not in the phylogeny
# as input, it takes model.pgls.traits_geomean (a list of one or more models of individual measurements on size, with names of the list elements matching the names of the columns), data.new (the data to size-scale, which can include data from species with and without phylogenetic information), scale_cols (vector of numbers of columns to size-scale), return_models (boolean to return the model fits), keep_geomean (boolean to retain the geomean column in the returned data frame)
# it returns either a data frame with the size-scaled data or a list containing the data frame with size-scaled data and the model fits
sizeScalePGLSNoPhylo <- function(model.pgls.traits_geomean, data.new, scale_cols, return_models=FALSE, keep_geomean=FALSE) {
  # check if models and data to scale are the same, and stop if they aren't
  if (!setequal(names(model.pgls.traits_geomean), colnames(data.new)[scale_cols])) {
    cat("Models and data to be scaled do not match.")
    return()
  }
  
  if (length(scale_cols)==1) {
    cat("Cannot scale a single column of data.")
    return()
  }
  
  scale_cols <- scale_cols[order(scale_cols)]  # put columns to scale in order
  model.pgls.traits_geomean <- model.pgls.traits_geomean[match(colnames(data.new)[scale_cols], names(model.pgls.traits_geomean))]  # sort the list of models to be in the same order as the columns to scale
  
  data.new$geomean <- apply(data.new[,scale_cols], MARGIN=1, mean)  # add a column containing the arithmetic mean of the log-transformed values (which is the same as the log of the geomean of the untransformed values)
  data.new.geomean_scaled <- data.new  # create a copy of the data frame to store size-scaled values

  for (i in 1:length(model.pgls.traits_geomean)) {  # loop over the models to use
    data.tmp.predicted <- predict(model.pgls.traits_geomean[[i]], newdata=data.frame(geomean=data.new.geomean_scaled$geomean))  # calculate the predicted values from the model
    data.new.geomean_scaled[,scale_cols[i]] <- data.new[,scale_cols[i]] - data.tmp.predicted  # calculate the residuals (the difference between the actual values and the predicted values), and store them in the data frame
  }
  
  if (!keep_geomean) data.new.geomean_scaled$geomean <- NULL
  if (return_models) {
    return(list(data.geomean_scaled=data.new.geomean_scaled, model.pgls.traits_geomean=model.pgls.traits_geomean))
  } else return(data.new.geomean_scaled)
}

## geomean-scale all the data, iterating over my list of variants
picidae.morph.log.fully_reduced.geomean_scaled.inc_no_phylo <- list()
for (i in names(picidae.morph.log.fully_reduced.geomean_scaled)) {  # loop over individual inclusion
  picidae.morph.log.fully_reduced.geomean_scaled.inc_no_phylo[[i]] <- sizeScalePGLSNoPhylo(model.pgls.traits_geomean = picidae.morph.log.fully_reduced.geomean_scaled[[i]]$model.pgls.traits_geomean, data.new=picidae.morph.log.fully_reduced[[i]]$data.imputed, scale_cols=2:ncol(picidae.morph.log.fully_reduced[[i]]$data.imputed))
}
rm(i)


## the function rotatePhylPCAnoPhylo() applies the rotation from phyl.pca data with or without phylogenetic information
# phyl.pca by default standardizes by the phylogenetic mean and does not standardize by the variances/sds, so this does, too
# as input, it takes pca.orig (the PCA generated previously), data.orig (the original data used to generate the PCA, which is necessary in order to properly apply the rotation), pca_cols.orig (the columns from data.orig that were included in the original PCA), phy.orig (the phylogeny used in the original PCA), method.orig (method passed to phyl.pca in the original PCA, necessary for calculating the phylogenetic mean), name_col.orig (name or number of the column in data.orig containing the taxon names), data.new (the data to apply the rotation to), pca_cols.new (the columns of data.new to apply the rotation to), name_col.new (name or number of the column in data.new containing the taxon names)
rotatePhylPCAnoPhylo <- function(pca.orig, data.orig, pca_cols.orig=NULL, phy.orig, method.orig="lambda", name_col.orig="My_genus_species", data.new, pca_cols.new=NULL, name_col.new="My_genus_species") {
  
  # convert data.orig and data.new to matrices with taxon names as rownames
  if (is.matrix(data.orig)) {
    data.orig.matrix <- data.orig
  } else if (is.data.frame(data.orig)) {
    data.orig.matrix <- data.matrix(data.orig[,pca_cols.orig])
    rownames(data.orig.matrix) <- data.orig[,name_col.orig]
  }
  if (is.matrix(data.new)) {
    data.new.matrix <- data.new
  } else if (is.data.frame(data.new)) {
    data.new.matrix <- data.matrix(data.new[,pca_cols.new])
    rownames(data.new.matrix) <- data.new[,name_col.new]
  }
  
  C <- vcv.phylo(phy.orig)[rownames(data.orig.matrix),rownames(data.orig.matrix)]  # get the phylogenetic variance-covariance matrix (necessary for calculating the phylogenetic mean)
  # calculate the phylogenetic mean of the original data
  if (method=="lambda") {
    mean.adjust <- phyl.vcv(X=data.orig.matrix, C=C, lambda=pca.orig$lambda)$alpha
  } else if (method=="BM") {
    mean.adjust <- phyl.vcv(X=data.orig.matrix, C=C, 1)$alpha
  }
  return(sweep(data.new.matrix, MARGIN=2, mean.adjust, FUN="-") %*% pca.orig$Evec)  # subtract the phylogenetic mean and multiply by the eigenvector matrix
}

## generate versions of unscaled data rotated by the phylogenetic pca, including data from species without phylogenetic information
picidae.morph.log.fully_reduced.phyl_pca.inc_no_phylo <- list()
for (i in names(picidae.morph.log.fully_reduced)) {  # loop over individual inclusion
  picidae.morph.log.fully_reduced.phyl_pca.inc_no_phylo[[i]] <- rotatePhylPCAnoPhylo(pca.orig = picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca, data.orig = picidae.morph.log.fully_reduced.treedata[[i]]$data, pca_cols.orig = 1:ncol(picidae.morph.log.fully_reduced.treedata[[i]]$data), phy.orig=picidae.morph.log.fully_reduced.treedata[[i]]$phy, data.new=picidae.morph.log.fully_reduced[[i]]$data.imputed, pca_cols.new=2:ncol(picidae.morph.log.fully_reduced[[i]]$data.imputed), method="lambda", name_col.new="My_genus_species")
}
rm(i)

## generate versions of size-scaled shape data rotated by the phylogenetic pca, including data from species without phylogenetic information
picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca.inc_no_phylo <- list()
for (i in names(picidae.morph.log.fully_reduced.geomean_scaled)) {  # loop over individual inclusion
  picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca.inc_no_phylo[[i]] <- rotatePhylPCAnoPhylo(pca.orig = picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca, data.orig = picidae.morph.log.fully_reduced.geomean_scaled[[i]]$data.geomean_scaled, pca_cols.orig = 2:ncol(picidae.morph.log.fully_reduced.geomean_scaled[[i]]$data.geomean_scaled), phy.orig=picidae.morph.log.fully_reduced.treedata[[i]]$phy, data.new=picidae.morph.log.fully_reduced.geomean_scaled.inc_no_phylo[[i]], pca_cols.new=2:ncol(picidae.morph.log.fully_reduced.geomean_scaled.inc_no_phylo[[i]]), method="lambda", name_col.orig="My_genus_species", name_col.new="My_genus_species")
}
rm(i)


### generate a matrix of pairwise Euclidean distances in pca space, and scale the pairwise overlap data by those distances (requires assuming morphology of taxa that don't have morphometric data (using phylogenetic or taxonomic proxies)

## identify proxies for taxa that don't have morphometric data (using phylogenetic or taxonomic proxies)
# note that I only use this proxy data to calculate similarity and similarity-scaled overlap aggregate scores for the taxa that do have data
# need different proxies for more inclusive/exclusive data sets; only bother with the ones I actually really plan to use (see freemind doc)
picidae.taxon.list.as_species <- read.csv(file="Picidae_taxon_list_as_species.csv",header=TRUE) # read in list of all taxa to generate complete matrix

picidae.taxon.list.proxies <- list()
picidae.taxon.list.proxies[["all_inds"]][["inc_no_phylo"]] <- picidae.taxon.list.as_species  # start with the full list of species
picidae.taxon.list.proxies[["all_inds"]][["inc_no_phylo"]]$proxy <- character(nrow(picidae.taxon.list.as_species))  # add empty column for proxy taxon
picidae.taxon.list.proxies[["all_inds"]][["inc_no_phylo"]]$proxy[picidae.taxon.list.proxies[["all_inds"]][["inc_no_phylo"]]$My_genus_species %in% rownames(picidae.morph.log.fully_reduced.phyl_pca.inc_no_phylo[["all_inds"]])] <- as.character(picidae.taxon.list.proxies[["all_inds"]][["inc_no_phylo"]]$My_genus_species[picidae.taxon.list.proxies[["all_inds"]][["inc_no_phylo"]]$My_genus_species %in% rownames(picidae.morph.log.fully_reduced.phyl_pca.inc_no_phylo[["all_inds"]])])  # if the taxon has data in the relevant object, set the proxy to be the taxon itself

picidae.taxon.list.proxies[["all_inds"]][["ex_no_phylo"]] <- picidae.taxon.list.as_species  # start with the full list of species
picidae.taxon.list.proxies[["all_inds"]][["ex_no_phylo"]]$proxy <- character(nrow(picidae.taxon.list.as_species))  # add empty column for proxy taxon
picidae.taxon.list.proxies[["all_inds"]][["ex_no_phylo"]]$proxy[picidae.taxon.list.proxies[["all_inds"]][["ex_no_phylo"]]$My_genus_species %in% rownames(picidae.morph.log.fully_reduced.phyl_pca[["all_inds"]]$pca$S)] <- as.character(picidae.taxon.list.proxies[["all_inds"]][["ex_no_phylo"]]$My_genus_species[picidae.taxon.list.proxies[["all_inds"]][["ex_no_phylo"]]$My_genus_species %in% rownames(picidae.morph.log.fully_reduced.phyl_pca[["all_inds"]]$pca$S)])  # if the taxon has data in the relevant object, set the proxy to be the taxon itself

for (i in names(picidae.taxon.list.proxies)) {  # loop over individual inclusion
  for (j in names(picidae.taxon.list.proxies[[i]])) {  # loop over include or exclude data without phylogenetic information
    write.csv(picidae.taxon.list.proxies[[i]][[j]], row.names=FALSE, file=paste("Picidae_taxon_list_as_species_proxies_", i, "_", j, ".csv", sep="")) # output the incomplete list (with proxies only for species with data) to a CSV; I then modified this in Excel
  }
}

table(picidae.morph[["raw.exclude_bad_juvs"]]$My_genus_species[grep("Veniliornis", picidae.morph[["raw.exclude_bad_juvs"]]$My_genus_species)])  # this allows screening the morph data to see how much coverage there is for a given species

## read in the proxy information after modifying the CSV file
for (i in names(picidae.taxon.list.proxies)) {  # loop over individual inclusion
  for (j in names(picidae.taxon.list.proxies[[i]])) {  # loop over include or exclude data without phylogenetic information
    picidae.taxon.list.proxies[[i]][[j]] <- read.csv(paste("Picidae_taxon_list_as_species_proxies_", i, "_", j, ".csv", sep=""), header=TRUE, stringsAsFactors=FALSE) # read in the complete list (with proxies for all species)
    picidae.taxon.list.proxies[[i]][[j]]$proxy_num <- match(picidae.taxon.list.proxies[[i]][[j]]$proxy, picidae.taxon.list.proxies[[i]][[j]]$My_genus_species) # add a column with the number of the proxy (in the list); this eases some downstream stuff
  }
}
rm(i,j)

## build a matrix of pairwise Euclidean distances in the unscaled pca space (as 1/(1+dist) between the taxa (or their proxies) in unscaled PCA space)
# if I don't use the inc_no_phylo data sets, I can only include taxa for which I have both genetic and morphological data
# I needed to make this bounded, or identical species (because of missing data) are going to make the values infinite.  I made it scale from 0 to 1, by calculating 1/(1+dist); that way if dist is infinite, the value is 0. if dist is 0, the value is 1

## the function calcSimilarityPairwise() generates a matrix of pairwise similarity (as 1/(1+dist)) for a set of traits, using a list of proxies
# as input, it takes data (as a matrix of trait values) and proxy_list (a data frame with columns for target taxon and proxy taxon)
calcSimilarityPairwise <- function(data, proxy_list) {
  if (is.vector(data)) data <- as.matrix(data)  # if data is a vector, convert it to a matrix, so that I can calculate distance for single traits
  pairwise.similarity <- matrix(nrow=nrow(proxy_list), ncol=nrow(proxy_list), dimnames=list(proxy_list$My_genus_species, proxy_list$My_genus_species))  # generate an empty matrix to fill in with the pairwise similarity values
  # loops through each taxon-by-taxon pair, and calculates pairwise similarity (as 1/(1+distance))
  for (m in 1:nrow(pairwise.similarity)) {  # loop over taxon 1
    for (n in 1:ncol(pairwise.similarity)) {  # loop over taxon 2
      pairwise.similarity[m,n] <- 1/(1+dist(rbind(data[proxy_list$proxy[m],], data[proxy_list$proxy[n],]))[1])  # rbinds the rows from the trait matrix for the current pair of taxa, and calculates 1/(1+distance) for that pair of taxa
    }
  }
  return(pairwise.similarity)
}

## calculate pairwise similarity matrices for geomean, unscaled_pca, and geomean_scaled_pca
picidae.morph.euclidean.pairwise <- list()
for (i in names(picidae.taxon.list.proxies)) {
  # for geomean
  picidae.morph.euclidean.pairwise[["geomean"]][[i]][["inc_no_phylo"]] <- calcSimilarityPairwise(data=picidae.morph.log.fully_reduced.geomean[[i]], proxy_list=picidae.taxon.list.proxies[[i]][["inc_no_phylo"]])
  
  # for unscaled_pca with inc_no_phylo
  picidae.morph.euclidean.pairwise[["phyl_pca"]][[i]][["inc_no_phylo"]] <- calcSimilarityPairwise(data=picidae.morph.log.fully_reduced.phyl_pca.inc_no_phylo[[i]], proxy_list=picidae.taxon.list.proxies[[i]][["inc_no_phylo"]])

  # for unscaled_pca with ex_no_phylo
  picidae.morph.euclidean.pairwise[["phyl_pca"]][[i]][["ex_no_phylo"]] <- calcSimilarityPairwise(data=picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S, proxy_list=picidae.taxon.list.proxies[[i]][["ex_no_phylo"]])

  # for geomean_scaled_pca with inc_no_phylo
  picidae.morph.euclidean.pairwise[["geomean_scaled.phyl_pca"]][[i]][["inc_no_phylo"]] <- calcSimilarityPairwise(data=picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca.inc_no_phylo[[i]], proxy_list=picidae.taxon.list.proxies[[i]][["inc_no_phylo"]])
  
  # for geomean_scaled_pca with ex_no_phylo
  picidae.morph.euclidean.pairwise[["geomean_scaled.phyl_pca"]][[i]][["ex_no_phylo"]] <- calcSimilarityPairwise(data=picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S, proxy_list=picidae.taxon.list.proxies[[i]][["ex_no_phylo"]])
}
rm(i)


### create an additional set of pairwise distribution overlap matrices, scaled by the Euclidean distances in geomean, unscaled_pca, and geomean_scaled_pca space
picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled <- list()
picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled <- list()
for (n in names(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0)) {  # loop over migratory vs. all_year distribution overlaps
  for (i in names(picidae.morph.euclidean.pairwise)) {  # loop over data variants
    for (j in names(picidae.morph.euclidean.pairwise[[h]][[i]])) {  # loop over individual inclusion
      for (k in names(picidae.morph.euclidean.pairwise[[h]][[i]][[j]])) {  # loop over inclusion/exclusion of data without phylogenetic information
        picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[n]][[i]][[j]][[k]] <- picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0[[n]][["overlaps.scaled"]] * picidae.morph.euclidean.pairwise[[i]][[j]][[k]] # multiply the two matrices element-wise (pairwise morphological similarity by the pairwise range overlap) to get similarity-scaled overlap for each taxon pair
        picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[n]][[i]][[j]][[k]] <- rowSums(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[n]][[i]][[j]][[k]]) - diag(picidae.overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[n]][[i]][[j]][[k]]) # for the overlap scaled by focal taxon range and similarity in PCA space, calculate the sum for each row (each taxon) and subtract the diagonal elements from the row sums (because all diagonal elements are 1 and inflate the values)
      }
    }
  }
}
rm(n,i,j,k)

n <- "migratory"
i <- "phyl_pca"
i <- "all_inds"
k <- "inc_no_phylo"

# output histogram of the range size (in km^2)
png(file='picidae.overlaps.BirdLife.hist.overlap.abs.png',height=400, width=480)
hist(log(diag(picidae.rangesize.shp.BirdLife.UnaryUnion.buffer0[[n]]/1000000)), main=NULL, xlab="log Range Size (km^2)",ylab="Count", col="gray")
dev.off()

# output histogram of the range overlap (overlap summed across focal taxon range)
png(file='picidae.overlaps.BirdLife.hist.overlap.sums.png',height=400, width=480)
hist(log(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[n]][["overlaps"]]/1000000),main=NULL,xlab="log Total Range Overlap (km^2)",ylab="Count", col="gray")
dev.off()

# output histogram of the range overlap (overlap summed across focal taxon range), scaled by focal taxon range size
png(file='picidae.overlaps.BirdLife.hist.overlap.scaled.sums.png',height=400, width=480)
hist(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[n]][["overlaps.scaled"]],main=NULL,xlab="Total Range Overlap scaled by Range Size",ylab="Count", col="gray")
dev.off()

## output histograms of the range overlap, scaled by focal taxon range size and similarity in geomean, unscaled PCA space, and geomean-scaled PCA space, using the most likely version of the morphology data
for (i in c("geomean", "phyl_pca", "geomean_scaled.phyl_pca")) {
  png(file=paste('picidae.overlaps.BirdLife.hist.overlap.scaled.euclidean_scaled.', h, '.sum.png', sep="") ,height=400, width=480)
  hist(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[n]][[i]][[j]][[k]], main=NULL, xlab="Range Overlap Scaled by Range Size and Morphological Similarity",ylab="Count", col="gray")
  dev.off()
}
rm(i,j,k,n)


### continue from here

###### run the analyses for Picinae ######

# list of non-Picinae taxa to exclude using grep "Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus"

#### load trees, overlaps, other data, and generate basic functions
#####

# drop non-Picinae taxa from trees
picinae.RAxML.all.BEAST_calibrated <- drop.tip(picidae.RAxML.all.BEAST_calibrated, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.BEAST_calibrated$tip.label, value=TRUE))
picinae.RAxML.all.BEAST_calibrated.no_calibs <- drop.tip(picidae.RAxML.all.BEAST_calibrated.no_calibs, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.BEAST_calibrated.no_calibs$tip.label, value=TRUE))
picinae.RAxML.all.BEAST_calibrated.with_proxies <- drop.tip(picidae.RAxML.all.BEAST_calibrated.with_proxies, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label, value=TRUE))
picinae.RAxML.all.r8s_calibrated.ladderized <- drop.tip(picidae.RAxML.all.r8s_calibrated.ladderized, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.r8s_calibrated.ladderized$tip.label, value=TRUE))
picinae.RAxML.all.r8s_calibrated.no_calibs.ladderized <- drop.tip(picidae.RAxML.all.r8s_calibrated.no_calibs.ladderized, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.r8s_calibrated.no_calibs.ladderized$tip.label, value=TRUE))
picinae.RAxML.all.rooted.ladderized <- drop.tip(picidae.RAxML.all.rooted.ladderized, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.rooted.ladderized$tip.label, value=TRUE))
picinae.RAxML.all.rooted.ladderized.chronos.lambda1 <- drop.tip(picidae.RAxML.all.rooted.ladderized.chronos.lambda1, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.RAxML.all.rooted.ladderized.chronos.lambda1$tip.label, value=TRUE))
picinae.mrbayes.all.ladderized <- drop.tip(picidae.mrbayes.all.ladderized, tip=grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.mrbayes.all.ladderized$tip.label, value=TRUE))

#####

### read in the morphological data from the database (caliper measurements)
picinae.morph <- list()
picinae.morph$raw <- picidae.morph$raw[grep("Jynx|Sasia|Verreauxia|Picumnus|Nesoctites|Hemicircus", picidae.morph$raw$My_genus_species, invert=TRUE),]


### to complete analyses for Picinae, use same code as above, replacing picidae with picinae for all variable and file names


### get specific information out of these processed morphological data

i <- "all_inds"

## check quantity of data for Picidae and Picinae
nrow(picidae.morph.log.fully_reduced[[c(i,"data.imputed")]])
nrow(picinae.morph.log.fully_reduced[[c(i,"data.imputed")]])
str(picidae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed")]])
str(picinae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed")]])
sum(picidae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed","My_genus_species")]] %in% picidae.morph.log.fully_reduced[[c(i,"data.imputed","My_genus_species")]])
sum(picinae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed","My_genus_species")]] %in% picinae.morph.log.fully_reduced[[c(i,"data.imputed","My_genus_species")]])

## determine how many taxa have measurement data
str(picinae.morph.log.fully_reduced, max.level=3)
length(unique(picidae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed","My_genus_species")]]))
length(unique(picidae.morph_combined.log.reduced_var_inds_sex.imputed[[c(i,"My_genus_species")]]))
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed[[c(i,"My_genus_species")]]
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.rm_taxa[[c(i,"My_genus_species")]]

## make a table of quantity of individuals of each sex measured for only the species that are in the fully-reduced data set, for Picidae and for Picinae
picidae.table.individuals.included.main_dataset <- table(picidae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed")]][picidae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed")]][["My_genus_species"]] %in% picidae.morph.log.fully_reduced[[c(i,"data.imputed","My_genus_species")]],c("My_genus_species","Sex")])
write.csv(picidae.table.individuals.included.main_dataset, file="picidae.measured_individuals.table.main_dataset.csv")
picinae.table.individuals.included.main_dataset <- table(picinae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed")]][picinae.morph_combined.log.reduced_var_inds.imputed[[c(i,"data.imputed")]][["My_genus_species"]] %in% picinae.morph.log.fully_reduced[[c(i,"data.imputed","My_genus_species")]],c("My_genus_species","Sex")])
write.csv(picinae.table.individuals.included.main_dataset, file="picinae.measured_individuals.table.main_dataset.csv")


## quantify the number of taxa (in Picidae and Picinae) for which I have all data (morph, in tree, and range)
picidae.morph.log.fully_reduced.treedata[[i]]$phy  # 134 species
picinae.morph.log.fully_reduced.treedata[[i]]$phy  # 121 species

## what effect did the inc_no_phylo have?
str(picidae.morph.log.fully_reduced.geomean_scaled[[c(i,"data.geomean_scaled")]], max.level=2)  # 134 species without inc_no_phylo
str(picinae.morph.log.fully_reduced.geomean_scaled[[c(i,"data.geomean_scaled")]], max.level=2)  # 121 species without inc_no_phylo
picidae.morph.log.fully_reduced.geomean_scaled[[c(i,"data.geomean_scaled","My_genus_species")]]

str(picidae.morph.log.fully_reduced.geomean_scaled.inc_no_phylo[[i]], max.level=2)  # 151 species with inc_no_phylo (same as the morphological data set)
str(picinae.morph.log.fully_reduced.geomean_scaled.inc_no_phylo[[i]], max.level=2)  # 134 species with inc_no_phylo (same as the morphological data set)


## checking how many taxa are in the trees
picidae.RAxML.all.BEAST_calibrated  # 178 species
picinae.RAxML.all.BEAST_calibrated  # 161 species


## getting info on the PCA axes (and plotting them!)
picidae.morph.log.fully_reduced.phyl_pca[[i]]$var_percents
picinae.morph.log.fully_reduced.phyl_pca[[i]]$var_percents

picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$L
picinae.morph.log.fully_reduced.phyl_pca[[i]]$pca$L

picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$var_percents
picinae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$var_percents

picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$L
picinae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$L

## checking correlations of similarity-scaled overlaps with raw (rangesize-scaled) overlaps

str(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0)
str(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled)

for(q in c("geomean","phyl_pca","geomean_scaled.phyl_pca")) {
  print(cor(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[c("mytax","migratory","overlaps.scaled")]], picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[c("mytax","migratory",q,i,j,k,l,"inc_no_phylo")]]))
}
rm(q)

pdf(file="/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Dissertation_chapters/Chapter_3_Diversification_Morph_Evol/Diss_ch3_Fig_???_overlaps_histograms.pdf")
par(mfrow=c(2,2))
hist(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[c("mytax","migratory","overlaps.scaled")]], col="gray", breaks=0:ceiling(max(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[c("mytax","migratory","overlaps.scaled")]])), xlab="Overlap scaled by focal taxon range size", main=NULL, xlim=c(0,20))
xlabs = c("Overlap scaled by average size", "Overlap scaled by unscaled morphologicla data", "Overlap scaled by size-scaled shape data")
names(xlabs) <- c("geomean","phyl_pca","geomean_scaled.phyl_pca")
for(q in c("geomean","phyl_pca","geomean_scaled.phyl_pca")) {
  hist(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[c("mytax","migratory",q,i,j,k,l,"inc_no_phylo")]], col="gray", breaks=0:ceiling(max(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[[c("mytax","migratory",q,i,j,k,l,"inc_no_phylo")]])), xlab=xlabs[q], main=NULL, xlim=c(0,20))
}
rm(q, xlabs)
par(opar)
dev.off()

# histogram of just the overlaps
pdf(file="/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Dissertation_chapters/Chapter_3_Diversification_Morph_Evol/Diss_ch3_Fig_2_overlaps_histograms_draft2.pdf")
hist(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[c("mytax","migratory","overlaps.scaled")]], col="gray", breaks=0:ceiling(max(picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[[c("mytax","migratory","overlaps.scaled")]])), xlab="Overlap scaled by focal taxon range size", main=NULL)
dev.off()

### plot the data by species, using standard scatterplots and phylomorphospace plots

## generate a data frame that differentiates Picidae and Picinae for use in plotting; used 1s (for Picinae) and 2s (for non-Picinae Picidae)
taxa.picinae_picidae <- read.csv(file="taxa.picinae_picidae.csv", row.names=1)  
names.tmp <- rownames(taxa.picinae_picidae)
taxa.picinae_picidae <- taxa.picinae_picidae[,1]
names(taxa.picinae_picidae) <- names.tmp
rm(names.tmp)

i <- "all_inds"; # start with the most basic version

## output scatterplot of PC1 and PC2 of unscaled data
pdf(file='picidae_unscaled_PC2_vs_PC1_plot.pdf', useDingbats=FALSE)
plot(picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S[,2] ~ picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S[,1], xlab=paste("PC1 (", format(round(picidae.morph.log.fully_reduced.phyl_pca[[i]]$var_percents[1],1), nsmall=1), "%)", sep=""), ylab=paste("PC2 (", format(round(picidae.morph.log.fully_reduced.phyl_pca[[i]]$var_percents[2],1), nsmall=1), "%)", sep=""), pch=19, col=c("black","gray60")[taxa.picinae_picidae[rownames(picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S)]]) # plot of PC1 and PC2 of unscaled morph data
legend(x="bottomleft", legend=c("Picinae", "non-Picinae"), pch=19, col=c("black","gray60"))
dev.off()

## output phylomorphospace of PC1 and PC2 of unscaled data
pdf(file='picidae_unscaled_PC2_vs_PC1_phylomorphospace.pdf', useDingbats=FALSE)
colors.tmp <- c(c("black","gray60")[taxa.picinae_picidae[rownames(picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S)]],rep("black",picidae.morph.log.fully_reduced.treedata[[i]]$phy$Nnode))
names(colors.tmp) <- 1:length(colors.tmp)
phylomorphospace(tree=picidae.morph.log.fully_reduced.treedata[[i]]$phy, X=picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S[,1:2], xlab=paste("PC1 (", format(round(picidae.morph.log.fully_reduced.phyl_pca[[i]]$var_percents[1],1), nsmall=1), "%)", sep=""), ylab=paste("PC2 (", format(round(picidae.morph.log.fully_reduced.phyl_pca[[i]]$var_percents[2],1), nsmall=1), "%)", sep=""), label="off", node.size=c(0,1), control=list(col.node=colors.tmp)) # phylomorphospace of PC1 and PC2 of geomean-scaled morph data
legend(x="bottomleft", legend=c("Picinae", "non-Picinae"), pch=19, col=c("black","gray60"))
rm(colors.tmp)
dev.off()

## output scatterplot of PC1 and PC2 of size-scaled shape data
pdf(file='picidae_geomean_scaled_PC2_vs_PC1_plot.pdf', useDingbats=FALSE)
plot(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S[,2] ~ picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S[,1], xlab=paste("PC1 (", round(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$var_percents[1],1), "%)", sep=""), ylab=paste("PC2 (", round(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$var_percents[2],1), "%)", sep=""), pch=19, col=c("black","gray60")[taxa.picinae_picidae[rownames(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S)]]) # phylomorphospace of PC1 and PC2 of geomean-scaled morph data
legend(x="bottomleft", legend=c("Picinae", "non-Picinae"), pch=19, col=c("black","gray60"))
dev.off()

## output phylomorphospace of PC1 and PC2 of size-scaled shape data
pdf(file='picidae_geomean_scaled_PC2_vs_PC1_phylomorphospace.pdf', useDingbats=FALSE)
colors.tmp <- c(c("black","gray60")[taxa.picinae_picidae[rownames(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S)]],rep("black",picidae.morph.log.fully_reduced.treedata[[i]]$phy$Nnode))
names(colors.tmp) <- 1:length(colors.tmp)
phylomorphospace(tree=picidae.morph.log.fully_reduced.treedata[[i]]$phy, X=picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S[,1:2], xlab=paste("PC1 (", round(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$var_percents[1],1), "%)", sep=""), ylab=paste("PC2 (", round(picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$var_percents[2],1), "%)", sep=""), label="off", node.size=c(0,1), control=list(col.node=colors.tmp)) # phylomorphospace of PC1 and PC2 of geomean-scaled morph data
legend(x="bottomleft", legend=c("Picinae", "non-Picinae"), pch=19, col=c("black","gray60"))
rm(colors.tmp)
dev.off()


### save off objects needed for other analyses

save(file='Picidae_morph_data_for_community_trait_analyses.RData', list=c("picidae.morph_combined.log.reduced_var_inds.imputed")) # save data needed for community trait distribution analyses

save(file='Picidae_data_for_distribution_morphology_evolution.RData', list=c(grep("overlaps|rangesize", ls(), value=TRUE), grep("all", ls(), value=TRUE), grep("fully_reduced", ls(), value=TRUE), grep("euclidean.pairwise", ls(), value=TRUE))) # save data needed for later analyses (to be used by Picidae_distribution_morphology_scripts)
