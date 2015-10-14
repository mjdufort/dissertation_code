## This is is one of several files containing scripts and functions used in processing and analysis of data for Matthew Dufort's Ph.D. dissertation at the University of Minnesota, titled "Coexistence, Ecomorphology, and Diversification in the Avian Family Picidae (Woodpeckers and Allies)."

## this file contains scripts for analyzing patterns of trait variation within and among communities of woodpeckers
## note: "geo cluster" and "community" are used interchangeably within this file; they have the same meaning


### load packages and data, and generate basic functions

## load necessary packages
library(caper)
library(ape)
library(plyr)
library(phytools)
library(geiger)
library(FNN)
library(fBasics)
library(nlme)
library(semTools)
library(gplots)
library(xlsx)
library(ggbiplot)
library(gplots)
library(lattice)

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

## pull in the morphology data from the Picidae_morphology_analyses R scripts
load(file='Picidae_morph_data_for_community_trait_analyses.RData')

## change "Colaptes_cafer" to "Colaptes_auratus" for all geo_clusters, with different data variations
for (i in 1:length(picidae.morph_combined.log.reduced_var_inds.imputed)) {
  picidae.morph_combined.log.reduced_var_inds.imputed[[i]]$My_genus_species <- gsub("Colaptes_cafer", "Colaptes_auratus", picidae.morph_combined.log.reduced_var_inds.imputed[[i]]$My_genus_species)
}
rm(i)

## pull in the phylogenetic data needed from Tree_manipulation.R
load(file='R_trees_for_combined_analyses.RData')

## import the geo_cluster membership data (from csv file)
picidae.geo_clusters.specimen_membership <- read.csv(file='Picidae_geo_cluster_membership_for_community_analyses.csv', header=TRUE, stringsAsFactors=FALSE)


### merge the morphological data and the community membership data

## generate translators to get full names of clusters and data variants from shorthand names
picidae.geo_clusters.to_use <- c("FL_Alachua","KS_Douglas","MI_Washtenaw","MN_Anoka","OR_JKL","WA_King")  # store vector of cluster shortnames
picidae.geo_clusters.name_translation <- c("Florida","Kansas","Michigan","Minnesota","Oregon","Washington")  # store vector of cluster full names, for using the names in figures and such  
names(picidae.geo_clusters.name_translation) <- picidae.geo_clusters.to_use  # set names of elements of full name vector to cluster shortnames
spec_inc.name_translation <- c("all focal populations", "all species", "all community members")  # store vector of full names of analysis variants
names(spec_inc.name_translation) <- c("all_focal_pops", "all_species", "comm_members")  # add shortnames as names of vector of full names of analysis variants

## the function mergeGeoClustersMorphData() merges the morphological data with the geo cluster membership data, with various defaults and filters
# as input, it takes morphdata (the morphological data, as a data frame), geo_clusters_specimen_membership (the cluster memberships, as a data frame), and variuos options to allow inclusion of specific geographic clusters, migratory vs. non-migratory species, ages, common vs. rare species, and extent of the cluster
# it returns a new data frame with a row for each matching pair of elements in the community membership data and morphological data
mergeGeoClustersMorphData <- function(morphdata, geo_clusters.specimen_membership, geo_clusters.to_use=c("FL_Alachua","KS_Douglas","MI_Washtenaw","MN_Anoka","OR_JKL","WA_King"), migratory_status.to_use=c("nm", "m_b"), ages.to_use="ad", species_status.to_use="c", locality_distance.to_use=c("in","ex")) {
  data.merged <- merge(x=subset(geo_clusters.specimen_membership, subset=((geo_cluster_shortname %in% geo_clusters.to_use) & (migratory_status %in% migratory_status.to_use) & (age_code %in% c(ages.to_use,"?")) & (species_status %in% species_status.to_use) & (locality_distance %in% locality_distance.to_use)), select=c("inst_catnum","geo_cluster_shortname")), y=morphdata, by.x="inst_catnum", by.y="Inst_CatNum")  # merge the data frames
  return(data.merged)
}

## generate variants of community morphological data
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters <- list()  # initialize a list to store variants of community morphological data
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds <-  merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership)  # merge the geo_cluster membership with the morphological data
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$complete_ind_only <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$complete_ind_only, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership)  # merge the geo_cluster membership with the morphological data, for complete individuals only
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_rare <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership, species_status.to_use=c("c","r"))  # merge the geo_cluster membership with the morphological data, with rare species included
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.ex_ex <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership, locality_distance.to_use="in")  # merge the geo_cluster membership with the morphological data, excluding individuals outside the geo_cluster circle
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$complete_ind_only.ex_ex <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$complete_ind_only, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership, locality_distance.to_use="in")  # merge the geo_cluster membership with the morphological data, for complete individuals only, excluding individuals outside the geo_cluster circle
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_m_nb <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership, migratory_status.to_use = c("nm", "m_b", "m_nb"))  # merge the geo_cluster membership with the morphological data, including probable non-breeders from migratory species
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_juv <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership, ages.to_use = c("ad","juv"))  # merge the geo_cluster membership with the morphological data, including juveniles
picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_juv_rare_m_nb <- merge_geo_clusters_morph_data(morphdata = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, geo_clusters.specimen_membership = picidae.geo_clusters.specimen_membership, ages.to_use = c("ad","juv"), migratory_status.to_use= c("nm", "m_b", "m_nb"), species_status.to_use=c("c","r"))  # merge the geo_cluster membership with the morphological data, including juveniles, rare species, and probable non-breeders from migratory species


### generate data objects by species, for later use

picidae.morph_combined.log.reduced_var_inds_sex.imputed <- list()  # initialize a list to store data by sex x species
picidae.morph_combined.log.reduced_var_inds_sex.imputed$all_inds <- reduce_to_mean_by_cols(picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, reduce_colnames=c("My_genus_species", "Sex"), mean_cols=4:ncol(picidae.morph_combined.log.reduced_var_inds.imputed$all_inds))  # reduce data set to single row for each sex x species combination

table.picidae.morph_combined.taxa_by_sex <- list()
table.picidae.morph_combined.taxa_by_sex$all_inds <- table(picidae.morph_combined.log.reduced_var_inds_sex.imputed$all_inds[,1:2]) # generate table of taxon (row) by sex (column), so that I can check which taxa have males and females

## determine taxa missing data for either males or females
picidae.morph_combined.taxa_by_mf <- list()
picidae.morph_combined.taxa_by_mf$no_female_male <- rownames(table.picidae.morph_combined.taxa_by_sex$all_inds[table.picidae.morph_combined.taxa_by_sex$all_inds[,"F"]==0 | table.picidae.morph_combined.taxa_by_sex$all_inds[,"M"]==0,]) # generate a vector of taxa that have either no males OR no females OR both
picidae.morph_combined.taxa_by_mf$female_male <- rownames(table.picidae.morph_combined.taxa_by_sex$all_inds[table.picidae.morph_combined.taxa_by_sex$all_inds[,"F"]==1 & table.picidae.morph_combined.taxa_by_sex$all_inds[,"M"]==1,]) # generate a vector of taxa that have both males AND females

## determine taxa to exclude
picidae.taxa.exclude <- list()
picidae.taxa.exclude$not_species <- c("Picidae_misidentified", grep(pattern="/|_x_", rownames(table.picidae.morph_combined.taxa_by_sex$all_inds), value=TRUE)) # generate a vector of ambiguous taxa to remove from the morphological data set; this finds hybrids and ambiguous taxa based on the "/" and "_x_" patterns

spp_for_geo_cluster_models <- read.csv("Picidae_additional_species_for_null_models.csv", header=TRUE)  # read in set of taxa for possible use in species-level models

## generate morphological data sets containing species means for taxa for possible use in species-level models
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models <- list()
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds <- rbind(reduce_to_mean_by_cols(data = subset(picidae.morph_combined.log.reduced_var_inds_sex.imputed$all_inds, subset=((Sex %in% c("F","M")) & (My_genus_species %in% picidae.morph_combined.taxa_by_mf$female_male) & (My_genus_species%in% spp_for_geo_cluster_models$My_genus_species) & !(My_genus_species %in% picidae.taxa.exclude$not_species))), reduce_colnames="My_genus_species", mean_cols = 3:ncol(picidae.morph_combined.log.reduced_var_inds_sex.imputed$all_inds)), reduce_to_mean_by_cols(data = subset(picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, subset=(My_genus_species %in% picidae.morph_combined.taxa_by_mf$no_female_male) & (My_genus_species%in% spp_for_geo_cluster_models$My_genus_species) & !(My_genus_species %in% picidae.taxa.exclude$not_species)), reduce_colnames="My_genus_species", mean_cols=4:ncol(picidae.morph_combined.log.reduced_var_inds.imputed$all_inds)))
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds <- picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds[order(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds$My_genus_species),] # sort by My_genus_species, to make things simpler
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds[picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds$My_genus_species=="Campephilus_principalis", "Pygostyle_width"] <- 2.738319  # this is the only data point that I needed to populate from the models used to impute additional data in my morphology analyses


### generate functions to calculate metrics, simulate communities, and generate expected distributions of metrics based on those simulations

## the function sample.vec() is a simple wrapper for sample() that allows more flexible sampling from a character vector
sample.vec <- function(x, ...) x[sample(length(x), ...)]

## the function sampleTaxa() samples from a list of taxa, while excluding specified tuples of taxa for which only one member can occur in a single sample
# as input, it takes taxa (a character vector containing the taxa to sample from), ntaxa (an integer, the number of taxa to be included in the sample), excludes.tuples (a list of zero or more character vectors, each containing a set of taxa that cannot occur together in a single sample), includes (a character vector containing taxa that must be included in the sample), and max.iter (the maximum number of attempts to make before stopping)
sampleTaxa <- function(taxa, ntaxa, excludes.tuples=list(c("Sphyrapicus_ruber","Sphyrapicus_varius","Sphyrapicus_nuchalis"),c("Colaptes_auratus","Colaptes_cafer","Colaptes_chrysoides"),c("Melanerpes_aurifrons","Melanerpes_carolinus","Melanerpes_uropygialis")), includes=NULL, max.iter=10000) {
  valid = FALSE
  iter = 0
  
  if (!is.null(includes)) {
    for (exclude in excludes.tuples) {
      if (length(intersect(includes, exclude)) > 1) {
        cat("Error: invalid sample. Inclusion set contains more than one member of an exclusion set.")
        return()
      } else if (length(intersect(includes, exclude)) == 1) {
        taxa <- taxa[!((taxa %in% exclude) & !(taxa %in% includes))]  # if any required taxon is in an exclusion set, remove all other taxa in that exclusion set from the list of taxa to sample from
      }
    }
  }
    
  while (!valid & (iter <= max.iter)) {  # repeat until a valid sample is found
    valid <- TRUE
    iter <- iter + 1
    
    if (length(union(includes,taxa)) == ntaxa) {
      taxa.sampled <- union(includes,taxa)  # if there is only one possible sample, return it
    } else taxa.sampled <- c(includes, sample(taxa[!(taxa %in% includes)], size=ntaxa-length(includes), replace=FALSE))  # generate a proposed sample of taxa
    
    for (i in 1:length(excludes.tuples)) {
      if (length(intersect(taxa.sampled, excludes.tuples[[i]])) > 1) valid <- FALSE  # check the current sample against the list of excluded sets of taxa, and reject the current sample if more than one member of any set is included
    }
    
    if (!valid & (iter>=max.iter)) {
      cat("Error: reached maximum numer of iterations without finding a valid sample.")
      return()
    }
  }
  
  return(taxa.sampled)
}

## the function simulate1.permutation() generates a simulated community by sampling from the set of available populations, with the constraint that no more than one population of a single taxon can be included in a sample
# as input, it takes ntaxa (the number of taxa to include in the sample) and traits.obs (a matrix containing the observed trait values, with taxon names as the rownames)
# it returns a matrix of trait values, with the taxon names as the rownames
simulate1.permutation <- function(ntaxa, traits.obs) {
  taxa <- unique(rownames(traits.obs))  # generate a vector of unique taxon names to sample
  traits <- matrix(data=NA, nrow=ntaxa, ncol=ncol(traits.obs))  # initialize and empty matrix to store permutation of trait data
  taxa.sampled <- sampleTaxa(taxa=taxa, ntaxa=ntaxa)  # sample the appropriate number of taxa from the vector of unique taxon names
  for (i in 1:ntaxa) {
    traits[i,] <- traits.obs[sample.vec(which(rownames(traits.obs) == taxa.sampled[i]), size=1),]  # loop over taxa to sample, storing trait data from a single population of that taxon in the permuted trait matrix
  }
  rownames(traits) <- taxa.sampled  # set the rownames of the permuted trait matrix to the taxon names
  return(traits)
}

## the function simulate1.normdist() generates a set of simulated trait values for a single community, using a normal distribution for trait values; it would be simple to expand this to include other distributions
# as input, it takes ntaxa (the number of taxa for the simulated community), ntraits (the number of traits), and distribution (a matrix with traits in rows, mean for each trait in column 1, and sd for each trait in column 2)
# it returns a matrix of trait values, with taxa in rows and traits in columns
simulate1.normdist <- function(ntaxa, ntraits, distribution) {
  if (ntraits != nrow(distribution)) {
    cat("Error: number of traits specified does not match number of traits in distribution matrix.")
    return()
  }
  traits <- matrix(NA, nrow=ntaxa, ncol=ntraits)
  for (i in 1:ntraits) traits[,i] <- rnorm(ntaxa, mean=distribution[i,1], sd=distribution[i,2])
  return(traits)
}

## the function simulate1.normdist() generates a set of simulated trait values for a single community, using a simulated BM or OU trait evolution on a phylogenetic tree (usually an estimated tree of the community members)
# if values are provided for alpha and theta, it uses OU; otherwise it uses BM
# as input, it takes phy (the phylogeny), ntraits (the number of traits), sigmasq (a vector containing the sigma^2 parameter values for the traits), and optionally alpha (a vector of alpha parameters for OU trait simulation), theta (a vector of theta parameters for OU trait simulation), and ntaxa (the number of taxa, if not using all taxon in phy
# it returns a matrix of simulated trait values, with taxa in rows and traits in columns
simulate1.phylogenetic <- function(phy, ntraits, sigmasq, alpha=NULL, theta=NULL, ntaxa=NULL) {
  require(phytools)
  
  # check the number of traits against the length of the sigma^2 vector
  if (ntraits != length(sigmasq)) {
    cat("Error: number of traits specified does not match number of traits in parameter vector.")
    return()
  }
  
  # if given a number of taxa to use, trim the tree to a sample of the taxa
  if (!is.null(ntaxa)) {
    taxa <- sampleTaxa(taxa=phy$tip.label, ntaxa=ntaxa)
    phy <- drop.tip(phy=phy, tip=phy$tip.label[!(phy$tip.label %in% taxa)])
  }
  
  traits <- matrix(nrow=length(phy$tip.label), ncol=ntraits)
  
  for (i in 1:ntraits) traits[,i] <- fastBM(tree=phy, sig2=sigmasq[i], alpha=alpha[i], theta=theta[i])  # loop over traits, simulating trait values using the trait-specific parameter values 
  rownames(traits) <- phy$tip.label  # set taxon names as rownames of the trait matrix
  return(traits)
}

## the function simulate.traits() is a wrapper for the various trait simulation methods; it simplifies the function for repeated simulations
# as input, it takes method ("BM", "OU", "permutation", or "distribution", the method to use in simulation) and various parameters to be passed to those functions
simulate.traits <- function(method, ntaxa=NULL, ntraits=NULL, traits.obs=NULL, distribution=NULL, phy=NULL, sigmasq=NULL, alpha=NULL, theta=NULL) {
  # method can be any of "BM", "OU", "permutation", or "distribution"
  traits <- switch(method,
                   BM = simulate1.phylogenetic(phy, ntraits, sigmasq, alpha, theta, ntaxa),
                   OU = simulate1.phylogenetic(phy, ntraits, sigmasq, alpha, theta, ntaxa),
                   distribution = simulate1.normdist(ntaxa, ntraits, distribution),
                   permutation = simulate1.permutation(ntaxa, traits.obs)
  )
  return(traits)
}

## the function simulateCommMetrics() generates distributions of community metrics based on permuted or simulated communities, with a number of options
# as input, it takes nsim (an integer, the number of simulations), method (the method to use for generating communities; can be "BM", "OU", "distribution", "permutation"), ntaxa (number of taxa for simulated communities), traits.obs (a matrix or data frame of empirical trait data, with taxon names as rownames or in the column name_col), name_col (if traits.obs is a data frame, the name or number of the column with the taxon names), trait_cols (if traits.obs is a data frame, the numbers of the columsn with trait data), various parameters to pass to simulation methods, metric (a character vector containing the metrics to be calculated), return.sims (a boolean to return the simulated communities), na_omit (a boolean to omit cases with NAs in the observed trait data), and print_every (output to the screen after every print_every iterations)
# it returns either a matrix of simulated community metrics, with metrics in columns and simulations in rows, or a list containing the matrix of simulated community metrics and a list containing a trait matrix for each simulated community
simulateCommMetrics <- function(nsim, method, ntaxa=NULL, traits.obs=NULL, name_col="My_genus_species", trait_cols=NULL, distribution=NULL, phy=NULL, sigmasq=NULL, alpha=NULL, theta=NULL, metric=c("var.mst", "var.mst.scaled", "rat.mst", "mean.mst", "min.mst", "mean.nndist", "sum.dist", "mean.dist", "var.dist", "cv.dist", "min.dist", "FRic", "FEve", "FDiv", "FDis", "kurtosis"), return.sims=FALSE, na_omit=TRUE, print_every=100) {

  # if traits.obs is a data frame, convert it to a matrix with taxon names as rownames
  if (is.data.frame(traits.obs)) {
    traits.obs.names <- traits.obs[,name_col]
    traits.obs <- data.matrix(traits.obs[,trait_cols])
    rownames(traits.obs) <- traits.obs.names
  }
  
  if (na_omit & !is.null(traits.obs)) traits.obs <- na.omit(traits.obs)  # remove NAs if specified
  if (return.sims) data.simulated <- vector("list", nsim)  # initialize list to store simulations, if specified
  
  # calculate number of traits
  if (method %in% c("BM","OU")) {
    ntraits <- length(sigmasq)
  } else if (method == "permutation") {
    ntraits <- ncol(traits.obs)
  } else if (method == "distribution") {
    ntraits <- nrow(distribution)
  }
  
  # if number of taxa not specified, calculate it
  if (is.null(ntaxa)) {
    if (method %in% c("BM","OU")) {
      ntaxa <- length(phy$tip.label)
    } else if ((method=="permutation") & (!is.null(traits.obs))) {
      ntaxa <- nrow(traits.obs)
    } else if (method=="distribution") {
      cat("Error: ntaxa must be specified for distribution simulation.")
      return()
    }
  }
  
  # if too few traits for FDiv, remove it from list of metrics
  if ((ntraits < 2) & ("FDiv" %in% metric)) {
    cat("FDiv cannot be calculated with fewer than 2 traits.\n")
    metric <- metric[metric != "FDiv"]
  }
  
  # if sample insufficient for kurtosis, remove it from list of metrics
  if ((ntraits >= ntaxa) & ("kurtosis" %in% metric)) {
    cat("Mardia's multivariate kurtosis cannot be calculated with more variables than samples.\n")
    metric <- metric[metric != "kurtosis"]
  }
  
  scores <- matrix(nrow=nsim, ncol=length(metric), dimnames=list(NULL, metric))  # initialize matrix to score community metrics
  
  for (i in 1:nsim) {
    if (i %% print_every == 0) cat("Simulation:", i, "\n")
    data.tmp <- simulate.traits(method, ntaxa=ntaxa, ntraits=ntraits, traits.obs=traits.obs, distribution=distribution, phy=phy, sigmasq=sigmasq, alpha=alpha, theta=theta)  # simulate data
    if (is.null(rownames(data.tmp))) rownames(data.tmp) <- letters[1:nrow(data.tmp)]  # add dummy rownames, for use in dbFD (which requires rownames)
    scores[i,] <- calcCommMetrics(traits.obs=data.tmp, metric=metric, na_omit=na_omit, quiet=TRUE)  # calculate metrics for current iteration
    if (return.sims) data.simulated[[i]] <- data.tmp  # if returning simulations, store current iteration in list
  }

  if (return.sims) {
    return(list(scores=scores, data.simulated=data.simulated))
  } else {
    return(scores)
  }
}

## the function calcCommMetrics() calculates community metrics from community trait data
# as input, it takes traits.obs (a matrix (with taxon names as rownames), or a data frame, with trait values for each community member in rows), name_col (if traits.obs is a data frame, the number or name of the column with taxon names), trait_cols (if traits.obs is a data frame, the numbers of the columns with trait data), metric (a character vector with the metrics to calculate), na_omit (a boolean to remove cases with NAs), and quiet (a boolean to output relevant information to the screen)
# it returns a vector of the community metrics, with elements of the vector named
calcCommMetrics <- function(traits.obs,  name_col="My_genus_species", trait_cols=NULL, metric=c("var.mst", "var.mst.scaled", "rat.mst", "mean.mst", "min.mst", "mean.nndist", "sum.dist", "mean.dist", "var.dist", "cv.dist", "min.dist", "FRic", "FEve", "FDiv", "FDis", "kurtosis"), na_omit=TRUE, quiet=FALSE) {
  require(vegan)
  require(FD)
  require(semTools)
  
  # if traits.obs is a data frame, convert it to a matrix with taxon names as rownames
  if (is.data.frame(traits.obs)) {
    traits.obs.names <- traits.obs[,name_col]
    traits.obs <- data.matrix(traits.obs[,trait_cols])
    rownames(traits.obs) <- traits.obs.names
  }
  
  if (na_omit) traits.obs <- na.omit(traits.obs)  # remove NAs if specified
  ntraits = ncol(traits.obs)  # calculate number of traits
  
  # if too few traits for FDiv, remove it from list of metrics
  if ((ntraits < 2) & ("FDiv" %in% metric)) {
    if (!quiet) cat("FDiv cannot be calculated with fewer than 2 traits.\n")
    metric <- metric[metric != "FDiv"]
  }
  
  # if sample insufficient for kurtosis, remove it from list of metrics
  if ((ntraits >= ntaxa) & ("kurtosis" %in% metric)) {
    if (!quiet) cat("Mardia's multivariate kurtosis cannot be calculated with more variables than samples.\n")
    metric <- metric[metric != "kurtosis"]
  }
  
  scores <- rep(NA, length(metric))  # intiialize vector to store metrics
  names(scores) <- metric  # add names of metrics to scores vector
  
  if (length(grep("mst", metric)) > 0) mst.tmp <- spantree(dist(traits.obs, method="euclidean"))  # find the minimum spanning tree if needed
  if ((length(grep("dist", metric)) > 0) | ("var.mst.scaled" %in% metric)) dist.tmp <- dist(traits.obs, method="euclidean")  # find the distance matrix if needed
  
  # calculate specified metrics
  if ("var.mst" %in% metric) scores["var.mst"] <- var(mst.tmp$dist)
  if ("var.mst.scaled" %in% metric) scores["var.mst.scaled"] <- var(mst.tmp$dist) / max(dist.tmp)
  if ("rat.mst" %in% metric) scores["rat.mst"] <- min(mst.tmp$dist) / max(mst.tmp$dist)
  if ("mean.mst" %in% metric) scores["mean.mst"] <- mean(mst.tmp$dist)
  if ("min.mst" %in% metric) scores["min.mst"] <- min(mst.tmp$dist)
  if ("mean.nndist" %in% metric) scores["mean.nndist"] <- mean.nndist(traits.obs)
  if ("sum.dist" %in% metric) scores["sum.dist"] <- sum(dist.tmp)
  if ("mean.dist" %in% metric) scores["mean.dist"] <- mean(dist.tmp)
  if ("var.dist" %in% metric) scores["var.dist"] <- var(dist.tmp)
  if ("min.dist" %in% metric) scores["min.dist"] <- min(dist.tmp)
  if ("cv.dist" %in% metric) scores["cv.dist"] <- sqrt(var(dist.tmp)) / mean(dist.tmp)
  if (length(intersect(metric, c("FRic", "FEve", "FDiv", "FDis"))) > 0) {
    dbFD.tmp <- dbFD(traits.obs, w.abun=FALSE, stand.x=FALSE, messages=FALSE)
    if ("FRic" %in% metric) scores["FRic"] <- dbFD.tmp$FRic
    if ("FEve" %in% metric) scores["FEve"] <- dbFD.tmp$FEve
    if ("FDiv" %in% metric) scores["FDiv"] <- dbFD.tmp$FDiv
    if ("FDis" %in% metric) scores["FDis"] <- dbFD.tmp$FDis
  }
  if ("kurtosis" %in% metric) {
    if (ntraits == 1) { scores["kurtosis"] <- kurtosis(traits.obs)[1]
    } else scores["kurtosis"] <- mardiaKurtosis(traits.obs)[1]
  }
  
  return(scores)
}


## the function mean.nndist() calculates the mean nearest neighbor distance
# as input, it takes a vector or matrix of trait values
# it returns an integer
mean.nndist <- function(traits.obs)  {
  if (is.null(nrow(traits.obs))) traits.obs <- cbind(traits.obs)  # if traits.obs is a vector, convert it to a single-column matrix
  neighbor.distances <- as.matrix(dist(traits.obs, method="euclidean"))  # generate distance matrix
  nndist <- rep(NA,nrow(traits.obs))  # initialize vector to store nearest neighbor distances
  for (i in 1:nrow(traits.obs))  nndist <- min(neighbor.distances[i,-i])  # loop over taxa, calculating the distance to the nearest neighbor in trait space
  return(mean(nndist))  # return the mean nearest neighbor distance
}

## the function calcDisplacementNNDist() determines nearest neighbor in trait space (within a community) to the taxon average value for each taxon in that community, and calculates displacement of focal populations towards or away from those nearest neighbors
# as input, it takes traits.obs (a matrix of trait observations with the taxon names as rownames), traits.tax_avg (a matrix of taxon average trait values with the taxon names as rownames), and tol (a distance below which a point is considered identical)
# returns a data frame containing the nearest neighbors (the taxon names), the trait value for that nearest neighbor, the distance to that nearest neighbor, and the displacement toward (positive) or away from (negative) that nearest neighbor
calcDisplacementNNDist <- function(traits.obs, traits.tax_avg, tol=0.0000000001) {
  taxa <- rownames(traits.obs)  # extract taxon names
  ntaxa <- length(taxa)  # get the number or species
  taxon.nn <- rep(NA, ntaxa)  # initialize a character vector to store the nearest neighbors
  traitvalue.nn <- matrix(NA, nrow=ntaxa, ncol=ncol(traits.obs))  # initialize a matrix to store the trait values of the nearest neighbors
  distance.nn_from_avg <- rep(NA, ntaxa)  # initialize a vector for the distance of the nearest neighbor from the taxon average of the focal species
  distance.nn_from_fp <- rep(NA, ntaxa)  # initialize a vector for the nearest neighbor from the focal population of the focal species
  distance.fp_from_avg <- rep(NA, ntaxa)  # intialize a vector for the distance of the focal population from the taxon average of the focal species (needed if I'm going to calculate the probability of moving towards nearest neighbor given the distance moved)
  prob.negative <- rep(NA, ntaxa)  # initialize a vector for the probability that the displacement will be negative (away from the nearest neighbor)
  
  dist.tax_avg.tmp <- apply(traits.tax_avg, MARGIN=1, function(x) {sqrt(colSums((t(traits.obs) - x)^2))})  # calculate distance matrix for the distances between the taxon averages (the columns) and the community populations (the rows)
  diag(dist.tax_avg.tmp) <- Inf  # set the diagonals to infinity so that I don't get a species as its own nearest neighbor
  rownames(dist.tax_avg.tmp) <- colnames(dist.tax_avg.tmp) <- taxa  # name the rows and columns of the distance matrix
  
  for (i in 1:ntaxa) {
    taxon.nn[i] <- taxa[which.min(dist.tax_avg.tmp[,i])]  # identify the community nearest neighbor to the taxon average of the focal species
    if (max(abs(traits.tax_avg[taxa[i],] - traits.obs[taxa[i],])) < tol) {
      distance.nn_from_avg[i] <- distance.nn_from_fp[i] <- distance.fp_from_avg[i] <- prob.negative[i] <- NA # set the distances (and therefore the displacement) to NA if I only have one population of that species; this eliminates meaningless 0s
    } else {
      traitvalue.nn[i,] <- traits.obs[taxon.nn[i],]  # get the trait values of the community nearest neighbor to the taxon average of the focal species
      distance.nn_from_avg[i] <- dist.tax_avg.tmp[taxon.nn[i], i]  # get the distance from the community nearest neighbor to the taxon average of the focal species
      distance.nn_from_fp[i] <- sqrt(sum((traits.obs[i,] - traitvalue.nn[i,])^2))  # get the distance from the community nearest neighbor to the focal population of the focal species
      distance.fp_from_avg[i] <- sqrt(sum((traits.obs[i,] - traits.tax_avg[taxa[i],])^2))  # get the distance of the focal population from the taxon average of the focal species
      if (distance.fp_from_avg[i] > (2 * distance.nn_from_avg[i])) {
        prob.negative[i] <- 1  # if the focal population is more than double the distance from the taxon average than the nearest neighbor is from the taxon average, it's impossible for the displacement to be positive
      } else {
        prob.negative[i] <- 0.5 + (0.5 * ((2 * asin(distance.fp_from_avg[i]/(2*distance.nn_from_avg[i]))) / pi))  # calculate the probability that the displacement is negative
      }
    }
  }
  displacement.towards_nn.from_avg <- distance.nn_from_avg - distance.nn_from_fp  # calculate the difference between the distances (positive values show displacement toward the nearest neighbor, negative values away)
  return(data.frame(taxon.fp=taxa, taxon.nn=taxon.nn, displacement=displacement.towards_nn.from_avg, distance.n_from_avg=distance.nn_from_avg, distance.fp_from_avg=distance.fp_from_avg, prob.negative=prob.negative))
}

## the function calcDisplacementGravity() calculates the difference in "gravity" (the inverses of the distances to all other community members, raised to the power exp, and summed) between species mean trait values and population trait values
# positive values indicate displacement toward nearby species/populations in trait space, negative values indicate displacement away
# as input, it takes traits.obs (a matrix of trait observations with the taxon names as rownames), traits.tax_avg (a matrix of taxon average trait values with the taxon names as rownames), tol (a distance below which a point is considered identical), and exp (the power to raise the inverse distances to before summing)
# it returns a data frame containing the taxon names and the change in gravity
calcDisplacementGravity <- function(traits.obs, traits.tax_avg, tol=0.0000000001, exp=1) {
  taxa <- rownames(traits.obs)  # store the taxon names
  ntaxa <- length(taxa)  # get the number or species
  
  dist.tax_avg.tmp <- apply(traits.tax_avg, MARGIN=1, function(x) {sqrt(colSums((t(traits.obs) - x)^2))})  # generate a matrix of the distances between the taxon averages (the columns) and the community populations (the rows)
  diag(dist.tax_avg.tmp) <- Inf  # set the diagonals to infinity so that I don't get a species as its own nearest neighbor
  rownames(dist.tax_avg.tmp) <- colnames(dist.tax_avg.tmp) <- taxa  # name the rows and columns of the distance matrix
  
  dist.focal_pop.tmp <- apply(traits.obs, MARGIN=1, function(x) {sqrt(colSums((t(traits.obs) - x)^2))})  # generate a matrix of the distances between the taxon averages (the columns) and the community populations (the rows)
  diag(dist.focal_pop.tmp) <- Inf  # set the diagonals to infinity so that I don't get a species as its own nearest neighbor
  rownames(dist.focal_pop.tmp) <- colnames(dist.focal_pop.tmp) <- taxa  # name the rows and columns of the distance matrix
  
  gravity.tax_avg <- gravity.focal_pop <- numeric(ntaxa)  # initialize vectors to store gravity values
  
  # loop over focal populations
  for (i in 1:ntaxa) {
    if (max(abs(traits.tax_avg[taxa[i],] - traits.obs[taxa[i],])) < tol) {
      gravity.tax_avg[i] <- gravity.focal_pop[i] <- NA # set the distances (and therefore the displacement) to NA if I only have one population of that species; this eliminates meaningless 0s
    } else {
      gravity.tax_avg[i] <- sum(1/(dist.tax_avg.tmp[,taxa[i]]^exp))  # calculate the inverses of the distances of the focal population's species mean to all other community members species mean trait values, raised to the power exp, and summed
      gravity.focal_pop[i] <- sum(1/(dist.focal_pop.tmp[,taxa[i]]^exp))  # calculate the inverses of the distances of the focal population to all other community members population trait values, raised to the power exp, and summed
    }  
  }
  gravity.change <- gravity.focal_pop - gravity.tax_avg  # calculate the difference in "gravity"
  return(data.frame(taxon.fp=taxa, displacement=gravity.change))
}

## the function simulate.displacement() generates expected distributions of displacement values
# it simulates communities by sampling one population at random for each species
# as input, it takes nsim (the number of simulations), taxa (the taxa to include in the simulated communities), traits.obs (a data frame containing population-level trait data), traits.tax_avg (a data frame containing species mean trait data), method ("nndist" or "gravity", the metric to use when calculating displacement values), and print_every (the frequency to output to the console)
simulate.displacement <- function(nsim, taxa, traits.obs, traits.tax_avg, method="nndist", print_every=100) {
  taxa <- taxa[order(taxa)]  # sort the taxon names
  traits.obs <- traits.obs[order(traits.obs$My_genus_species),]  # sort the population-level trait data by taxon name
  traits.tax_avg <- traits.tax_avg[order(traits.obs$My_genus_species),]  # sort the species mean trait data by taxon name
  traits.tax_avg.tmp <- data.matrix(traits.tax_avg[traits.tax_avg$My_genus_species %in% taxa, 2:ncol(traits.tax_avg)])  # create a temporary species mean trait matrix including only the taxa to be included in simulated communities
  rownames(traits.tax_avg.tmp) <- traits.tax_avg$My_genus_species[traits.tax_avg$My_genus_species %in% taxa]  # set the rownames of the temporary species mean trait matrix to be the species names
  traits.obs.matrix <- data.matrix(traits.obs[,4:ncol(traits.obs)])  # extract the population-level trait data as a matrix
  rownames(traits.obs.matrix) <- traits.obs$My_genus_species  # set the rownames of the population-level trait matrix to be the species names
  displacements.simulated <- list()  # initialize a list to store the displacement values from the simulated communities
  
  for (i in 1:nsim) {
    if (i %% print_every == 0) cat("Simulation:", i, "\n")
    
    # create a simulated trait matrix by sampling one population from each species in the vector taxa
    traits.obs.tmp <- matrix(data=NA, nrow=length(taxa), ncol=(ncol(traits.obs)-3))
    for (j in 1:length(taxa)) {
      traits.obs.tmp[j,] <- traits.obs.matrix[sample.vec(which(rownames(traits.obs.matrix) == taxa[j]), size=1),]
    }
    rownames(traits.obs.tmp) <- taxa
    
    # calculate the displacement values for the simulated community, and store them in the list
    if (method=="nndist") {
      displacements.simulated[[i]] <- calcDisplacementNNDist(traits.obs=traits.obs.tmp, traits.tax_avg=traits.tax_avg.tmp)
    } else if (method=="gravity") {
      displacements.simulated[[i]] <- calcDisplacementGravity(traits.obs=traits.obs.tmp, traits.tax_avg=traits.tax_avg.tmp)
    }
  }
  
  return(displacements.simulated)
}


### generate data frames with the populations to be analyzed (unique combinations of species and geo_cluster)

species.geo_cluster.combinations <- list()  # initialize list to store variants of data

## generate a data frame with two columns, containing the unique combinations of species and geo_clusters for which I have data
species.geo_cluster.combinations$base <- unique(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds[c("My_genus_species","geo_cluster_shortname")])  # extract the unique combinations of species and geo_cluster
species.geo_cluster.combinations$base <- species.geo_cluster.combinations$base[order(species.geo_cluster.combinations$base$My_genus_species, species.geo_cluster.combinations$base$geo_cluster_shortname),]  # sort by species and geo_cluster name
species.geo_cluster.combinations$base$full_name <- apply(species.geo_cluster.combinations$base, MARGIN=1, FUN = function(x) paste(x[1],x[2], sep = "_"))  # add a column with the concatenated name

# generate a data frame with two columns, containing the unique combinations of species and geo_clusters for which I have data, including rare species, juveniles, and possible out-of-season individuals from migratory species
species.geo_cluster.combinations$inc_juv_rare_m_nb <- unique(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_juv_rare_m_nb[c("My_genus_species","geo_cluster_shortname")])  # extract the unique combinations of species and geo_cluster
species.geo_cluster.combinations$inc_juv_rare_m_nb <- species.geo_cluster.combinations$inc_juv_rare_m_nb[order(species.geo_cluster.combinations$inc_juv_rare_m_nb$My_genus_species, species.geo_cluster.combinations$inc_juv_rare_m_nb$geo_cluster_shortname),]  # sort by species and geo_cluster name
species.geo_cluster.combinations$inc_juv_rare_m_nb$full_name <- apply(species.geo_cluster.combinations$inc_juv_rare_m_nb, MARGIN=1, FUN = function(x) paste(x[1],x[2], sep="_"))  # add a column with the concatenated name


### check the data for outliers by plotting histograms of all individuals from each population

## output histograms of each measurement for each population
for (i in 1:nrow(species.geo_cluster.combinations$base)) {  # loop over populations
  pdf(file=paste(species.geo_cluster.combinations$base$My_genus_species[i],"_", species.geo_cluster.combinations$base$geo_cluster_shortname[i], "_histogram_matrix.pdf", sep=""), height=8, width=6)
  par(mfrow=c(4,3))
  for (j in 5:ncol(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds)) {  # loop over traits
    data.tmp <- subset(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds, subset=((My_genus_species==species.geo_cluster.combinations$base$My_genus_species[i]) & (geo_cluster_shortname==species.geo_cluster.combinations$base$geo_cluster_shortname[i])))[,j]  # extract the data for current population and trait
    if (any(!is.na(data.tmp))) {hist(data.tmp, main=colnames(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds)[j], breaks=10)  # output histogram if there is data
    } else plot.new()
  }
  dev.off()
  rm(data.tmp)
}
rm(i,j)


## output histograms of each measurement for each population, including rare species, juveniles, and possible out-of-season individuals from migratory species
for (i in 1:nrow(species.geo_cluster.combinations$inc_juv_rare_m_nb)) {  # loop over populations
  pdf(file=paste(species.geo_cluster.combinations$inc_juv_rare_m_nb$My_genus_species[i],"_", species.geo_cluster.combinations$inc_juv_rare_m_nb$geo_cluster_shortname[i], "_inc_rare_juv_m_nb_histogram_matrix.pdf", sep=""), height=8, width=6)
  par(mfrow=c(4,3))
  for (j in 5:ncol(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_juv_rare_m_nb)) {  # loop over measurements
    data.tmp <- subset(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_juv_rare_m_nb, subset=((My_genus_species==species.geo_cluster.combinations$inc_juv_rare_m_nb$My_genus_species[i]) & (geo_cluster_shortname==species.geo_cluster.combinations$inc_juv_rare_m_nb$geo_cluster_shortname[i])))[,j]  # extract data for current population and measurement
    if (any(!is.na(data.tmp))) {hist(data.tmp, main=colnames(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds.inc_juv_rare_m_nb)[j], breaks=10)  # if there is data, output a histogram
    } else plot.new()
  }
  dev.off()
  rm(data.tmp)
}
rm(i,j)


### take averages of all measurements by sex, then by species for each population

## find the mean of all measurements for each sex of each taxon in each geo_cluster, averaging across individuals
picidae.morph_combined.log.reduced_var_inds_sex.imputed.geo_clusters <- list()  # initialize list to store data by sex
picidae.morph_combined.log.reduced_var_inds_sex.imputed.geo_clusters$all_inds <-  reduce_to_mean_by_cols(data = picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds, reduce_colnames = c("geo_cluster_shortname", "My_genus_species", "Sex"), mean_cols = 5:ncol(picidae.morph_combined.log.reduced_var_inds.imputed.geo_clusters$all_inds))

## find the mean of all measurements for each taxon in each geo_cluster, averaging across sexes
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.geo_clusters <- list()  # initialize list to store data by species
picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.geo_clusters$all_inds <- reduce_to_mean_by_cols(data = subset(picidae.morph_combined.log.reduced_var_inds_sex.imputed.geo_clusters$all_inds, subset=(Sex %in% c("F","M"))), reduce_colnames = c("geo_cluster_shortname", "My_genus_species"), mean_cols = 4:ncol(picidae.morph_combined.log.reduced_var_inds_sex.imputed.geo_clusters$all_inds))

## create a fully reduced version of the morphological data
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters <- list()
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds <- na.omit(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.geo_clusters$all_inds)  # drop geo_cluster x species combinations where I don't have all data
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds$full_name <- apply(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds, MARGIN=1, FUN = function(x) paste(x[2],x[1], sep = "_")) # add a column to my collection of species with the concatenated name
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds <- picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds[c(1:2, ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds), 3:(ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds)-1))] # re-order the columns


### generate tree for use in population-level PGLS and phylogenetic PCA, by adding polytomies to the tree that split each species into the populations I have
## original species are not dropped, so that the same tree can be used later for populations or species

## add populations as a star tree of n taxa with edge.length set to something (0.2 for now) and root.edge set to the taxon's terminal edge length minus the edge length in the star tree
# this was done on the full tree, so that there is a version of the tree that can be trimmed to include any combination of population and species
picidae.RAxML.all.BEAST_calibrated.geo_clusters <- picidae.RAxML.all.BEAST_calibrated  # make a copy of the phylogeny, to which the polytomies will be added
picidae.RAxML.all.BEAST_calibrated.geo_clusters[grep("rate", names(picidae.RAxML.all.BEAST_calibrated.geo_clusters))] <- NULL # remove the extra variables that will nto be used
picidae.species.geo_clusters <- unique(species.geo_cluster.combinations$base$My_genus_species)

## find the depth of the shallowest node among species in my geo_clusters in the BEAST tree; this gives some idea of when the population polytomy should be dated
min(max(node.depth.edgelength(picidae.RAxML.all.BEAST_calibrated)) - node.depth.edgelength(picidae.RAxML.all.BEAST_calibrated)[picidae.RAxML.all.BEAST_calibrated$edge[picidae.RAxML.all.BEAST_calibrated$edge[,2] %in% which(picidae.RAxML.all.BEAST_calibrated$tip.label %in% picidae.species.geo_clusters),1]])

## identify the edges that lead to the tips that are species in these community analyses; also identify the nodes that are directly ascendant from those tips
which(picidae.RAxML.all.BEAST_calibrated$tip.label %in% picidae.species.geo_clusters)  # returns the node numbers of the species in community analyses
picidae.RAxML.all.BEAST_calibrated$edge[picidae.RAxML.all.BEAST_calibrated$edge[,2] %in% which(picidae.RAxML.all.BEAST_calibrated$tip.label %in% picidae.species.geo_clusters),1]  # returns the numbers of the ascendant nodes from the species in community analyses
which(picidae.RAxML.all.BEAST_calibrated$edge[,2] %in% which(picidae.RAxML.all.BEAST_calibrated$tip.label %in% picidae.species.geo_clusters))  # returns the numbers of the ascendant edges from the species in community analyses

## generate a list of star trees for the species in my geo_clusters, with branch lengths of 0.2, and root edge set to the length of the edge leading to that species in the original tree minus the length of the terminal edges in the star tree
picidae.star_trees.species.geo_clusters <- list() # create an empty list to store the star trees
picidae.star_tree.brlen <- 0.2  # set the length of the terminal edges in the star tree (in millions of years)
for (i in 1:length(picidae.species.geo_clusters)) {
  picidae.star_trees.species.geo_clusters[[i]] <- stree(n = sum(species.geo_cluster.combinations$base$My_genus_species == picidae.species.geo_clusters[i]), type="star", tip.label = species.geo_cluster.combinations$base$full_name[species.geo_cluster.combinations$base$My_genus_species == picidae.species.geo_clusters[i]])
  # picidae.star_trees.species.geo_clusters[[i]]$edge.length <- rep(max(picidae.star_tree.brlen, min((picidae.RAxML.all.BEAST_calibrated$edge.length[which(picidae.RAxML.all.BEAST_calibrated$edge[,2] == which(picidae.RAxML.all.BEAST_calibrated$tip.label == picidae.species.geo_clusters[i]))] / 4), 1.0)), nrow(picidae.star_trees.species.geo_clusters[[i]]$edge))  # variant to set edge.length = max(0.2 (=200,000 years) or 1/4 of the length of the ascendant edge)
  picidae.star_trees.species.geo_clusters[[i]]$edge.length <- rep(picidae.star_tree.brlen, nrow(picidae.star_trees.species.geo_clusters[[i]]$edge)) # set edge.length = 0.2 (=200,000 years)
  picidae.star_trees.species.geo_clusters[[i]]$root.edge <- picidae.RAxML.all.BEAST_calibrated$edge.length[which(picidae.RAxML.all.BEAST_calibrated$edge[,2] == which(picidae.RAxML.all.BEAST_calibrated$tip.label == picidae.species.geo_clusters[i]))] - picidae.star_trees.species.geo_clusters[[i]]$edge.length[1]  # add a root edge that is the length of the terminal edge in the original tree minus the length of the terminal edges in the star tree
  names(picidae.star_trees.species.geo_clusters)[i] <- picidae.species.geo_clusters[i]  # set the name of the current list element to the species
}
rm(i)

## bind the population-level star trees to the full tree
for (i in 1:length(picidae.star_trees.species.geo_clusters)) {
  picidae.RAxML.all.BEAST_calibrated.geo_clusters <- bind.tree(picidae.RAxML.all.BEAST_calibrated.geo_clusters, picidae.star_trees.species.geo_clusters[[i]], where=picidae.RAxML.all.BEAST_calibrated.geo_clusters$edge[picidae.RAxML.all.BEAST_calibrated.geo_clusters$edge[,2] == which(picidae.RAxML.all.BEAST_calibrated.geo_clusters$tip.label == names(picidae.star_trees.species.geo_clusters)[i]),1])
}
rm(i)
picidae.RAxML.all.BEAST_calibrated.geo_clusters <- collapse.singles(picidae.RAxML.all.BEAST_calibrated.geo_clusters)  # collapse single nodes, which can be introduced by the previous steps

## output the tree with population-level polytomies
pdf(file="RAxML_all_tree.BEAST_calibrated.rooted.ladderized.geo_clusters.pdf", width=4.5, height=10)
plot(picidae.RAxML.all.BEAST_calibrated.geo_clusters, direction="rightwards", cex=0.29, no.margin=T, x.lim=c(-7, 72), label.offset=0.2)
dev.off()

## create a copy of the tree with only the geo_clusters
picidae.RAxML.all.BEAST_calibrated.geo_clusters_only <- drop.tip(picidae.RAxML.all.BEAST_calibrated.geo_clusters, tip=picidae.RAxML.all.BEAST_calibrated.geo_clusters$tip.label[!(picidae.RAxML.all.BEAST_calibrated.geo_clusters$tip.label %in% species.geo_cluster.combinations$base$full_name)])  

## create a copy of the tree that only has the species included in the community analyses
picidae.RAxML.all.BEAST_calibrated.geo_clusters_species_only <- drop.tip(picidae.RAxML.all.BEAST_calibrated.geo_clusters, tip=picidae.RAxML.all.BEAST_calibrated.geo_clusters$tip.label[!(picidae.RAxML.all.BEAST_calibrated.geo_clusters$tip.label %in% picidae.species.geo_clusters)])


### scale data by the geomean (in various ways)

## the function sizeScale.pgls() takes a data frame of log-transformed data values, with values to scale in scale_cols, and a tree containing the taxa those data represent, and adds a geomean column to the end, then scales the variables in scale_cols by the geomean (by finding the residuals from PGLS regression of each variable on the geomean)
#  it can do this with incomplete overlap in taxa; it uses treedata to drop the taxa not found in both, and outputs a list of the taxa that were dropped
# as input, it takes data (a data frame containing species means), scale_cols (vector of numbers of columns to scale), phy (phylogeny), name_col (name or number of the column of the data frame with the taxon names), keep_geomean (boolean to retain the geomean column in the returned data frame), return_models (boolean to return the model fits), and optional methods to pass to the gls function for PGLS model fitting
# it returns either a data frame with the size-scaled data (with or without geomean) or a list containing the size-scaled data and the model fits
sizeScale.pgls <- function(data, scale_cols, phy, name_col, keep_geomean=FALSE, return_models=FALSE, method="ML", control=list()) {
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

## the function sizeScale.pgls.noPhylo() uses the models built from PGLS regression of measurements on size to scale data not included in those models because the species are not in the phylogeny
# as input, it takes model.pgls.traits_geomean (a list of one or more models of individual measurements on size, with names of the list elements matching the names of the columns), data.new (the data to size-scale, which can include data from species with and without phylogenetic information), scale_cols (vector of numbers of columns to size-scale), return_models (boolean to return the model fits), keep_geomean (boolean to retain the geomean column in the returned data frame)
# it returns either a data frame with the size-scaled data or a list containing the data frame with size-scaled data and the model fits
sizeScale.pgls.noPhylo <- function(model.pgls.traits_geomean, data.new, scale_cols, return_models=FALSE, keep_geomean=FALSE) {
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

## get the geomean as a separate object
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean <- list()
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds <- data.frame(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds[1:3], geomean=apply(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds[,4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds)],MARGIN=1,mean))

## create a list object to store all my geomean_scaled data
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled <- list()
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls_fits <- sizeScale.pgls(data = picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds, scale_cols=4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds), phy=picidae.RAxML.all.BEAST_calibrated.geo_clusters, name_col=3, return_models=TRUE)  # geomean-scale the data and return the models
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls <- picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls_fits$data.geomean_scaled  # extract just the data


### conduct phylogenetic PCA on populations, and calculate species means for the data variants

## generate a new data frame containing the species means (as the mean of all populations of that species)
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops <- list()
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops$all_inds <- reduce_to_mean_by_cols(data = picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds, reduce_colnames = "My_genus_species", mean_cols=4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds))



## the function wrapped.phyl_pca() is a wrapper to run phylogenetic PCA; it creates the necessary intermediate steps (matrix, treedata), and also calculates the percentage of the variance captured by each PCA
# as input, it takes data (the data as a data frame), pca_cols (the columns on which to conduct phylogenetic PCA), phy (the phylogeny), name_col (the name or number of the column containing taxon names), and method (the method to pass to phyl.pca)
# it returns a list containing the results of the phylogenetic PCA (an object of class phyl.pca) and a vector of the percent variance explained by each PCA axis
wrapped.phyl_pca <- function(data, pca_cols, phy, name_col, method="lambda") {
  require(geiger)
  require(phytools)
  newdata.matrix <- data.matrix(data[,pca_cols])  # create a matrix version of the data (in order to use geiger::treedata)
  rownames(newdata.matrix) <- data[,name_col]  # set the rownames to be the taxon names (treedata needs this)
  newdata.treedata <- treedata(phy = phy, data = newdata.matrix, warnings=!quiet, sort=TRUE)  # construct a treedata object (which reduces the phylogeny and the data matrix based on shared membership)
  newdata.pca <- phyl.pca(tree = newdata.treedata$phy, Y = newdata.treedata$data, method=method, mode="cov")  # generate a phylogenetic PCA of the values
  newdata.variances <- diag(newdata.pca$Eval)/sum(newdata.pca$Eval)*100  # calculate the percentage of the variance in the data that is captured be each PCA axis
  return(list(pca = newdata.pca, var_percents = newdata.variances))
}

## generate a phylogenetic PCA of the geo_cluster means, for the unscaled morphological data and the geomean-scaled data
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca <- list()
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds <- wrapped.phyl_pca(data = picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds, phy = picidae.RAxML.all.BEAST_calibrated.geo_clusters, pca_cols = 4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds), name_col = "full_name")
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds <- wrapped.phyl_pca(data = picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls, phy = picidae.RAxML.all.BEAST_calibrated.geo_clusters, pca_cols = 4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls), name_col = "full_name")

## the function rotate.phyl_pca.noPhylo() applies the rotation from phyl.pca to data with or without phylogenetic information
# phyl.pca by default standardizes by the phylogenetic mean and does not standardize by the variances/sds, so this does, too
# as input, it takes pca.orig (the PCA generated previously), data.orig (the original data used to generate the PCA, which is necessary in order to properly apply the rotation), pca_cols.orig (the columns from data.orig that were included in the original PCA), phy.orig (the phylogeny used in the original PCA), method.orig (method passed to phyl.pca in the original PCA, necessary for calculating the phylogenetic mean), name_col.orig (name or number of the column in data.orig containing the taxon names), data.new (the data to apply the rotation to), pca_cols.new (the columns of data.new to apply the rotation to), name_col.new (name or number of the column in data.new containing the taxon names)
rotate.phyl_pca.noPhylo <- function(pca.orig, data.orig, pca_cols.orig=NULL, phy.orig, method.orig="lambda", name_col.orig="My_genus_species", data.new, pca_cols.new=NULL, name_col.new="My_genus_species") {
  
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


### build up objects with species means, to use for null model means
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models <- list()
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds <- na.omit(picidae.morph_combined.log.reduced_var_inds_sex_taxa.imputed.spp_for_geo_cluster_models$all_inds)

## calculate the geomean for the species mean data
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.geomean <- list()
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.geomean$all_inds <- data.frame(picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds[1], geomean=apply(picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds[,2:ncol(picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds)],MARGIN=1,mean))

## geomean-scale the data for the spp_for_geo_cluster_models
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.geomean_scaled <- list()
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.geomean_scaled$all_inds.pgls <- sizeScale.pgls.noPhylo(model.pgls.traits_geomean = picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls_fits$model.pgls.traits_geomean, data.new = picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds, scale_cols=2:ncol(picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds))

## rotate the unscaled and the geomean_scaled.pgls data for the species means, by the PCAs of the populations, so that they use the same rotation
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.phyl_pca_rotated <- list()
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.phyl_pca_rotated$unscaled.all_inds <- rotate.phyl_pca(pca.orig=picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca, data.orig=picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds, pca_cols.orig=4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters$all_inds), phy.orig=picidae.RAxML.all.BEAST_calibrated.geo_clusters, data.new=picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds, pca_cols.new=2:ncol(picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models$all_inds), method="lambda", name_col.orig="full_name", name_col.new="My_genus_species")
picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.phyl_pca_rotated$geomean_scaled.pgls.all_inds <- rotate.phyl_pca(pca.orig=picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca, data.orig=picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls, pca_cols.orig=4:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls), phy.orig=picidae.RAxML.all.BEAST_calibrated.geo_clusters, data.new=picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.geomean_scaled$all_inds.pgls, pca_cols.new=2:ncol(picidae.morph_combined.log.fully_reduced.imputed.spp_for_geo_cluster_models.geomean_scaled$all_inds.pgls), method="lambda", name_col.orig="full_name", name_col.new="My_genus_species")


### various plots of geomean and PCA axes for the populations

## store values to use as the boundaries of the plots (so that they're consistent, even with different data)
xylims <- list()
xylims$unscaled.pc12.x <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,1])
xylims$unscaled.pc12.y <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,2])
# xylims$pgls_scaled.pc12.x <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1])
xylims$pgls_scaled.pc12.x <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1]) + c(0,0.05) # adjusted to allow space for the legend
xylims$pgls_scaled.pc12.y <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,2])
xylims$geomean.pgls_scaled.pc1.x <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$geomean)
xylims$geomean.pgls_scaled.pc1.y <- range(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1])

## plot PC1 vs. PC2 for each of the geo_clusters, using the phylogenetic PCA of the unscaled morphological data
for (i in picidae.geo_clusters.to_use) {
  pdf(file=paste("Picidae_community_unscaled_PC1_vs_PC2_", i, ".pdf", sep=""))
  with(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds, plot(pca$S[grep(i, rownames(pca$S)),1], pca$S[grep(i, rownames(pca$S)),2], xlab = paste("PC1 (", round(var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(var_percents[2], 2), "%)", sep=""), xlim=xylims$unscaled.pc12.x, ylim=xylims$unscaled.pc12.y))
  dev.off()
}

## plot PC1 vs. PC2 for each of the geo_clusters, using the phylogenetic PCA of the morphological data scaled by the geomean using residuals from PGLS regression of each variable on the geomean
for (i in picidae.geo_clusters.to_use) {
  pdf(file=paste("Picidae_community_geomean_scaled_pgls_PC1_vs_PC2_", i, ".pdf", sep=""), width=6, height=6)
  with(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds, plot(pca$S[grep(i, rownames(pca$S)),1], pca$S[grep(i, rownames(pca$S)),2], xlab = paste("PC1 (", round(var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(var_percents[2], 2), "%)", sep=""), xlim=xylims$pgls_scaled.pc12.x, ylim=xylims$pgls_scaled.pc12.y, pch=16))
  dev.off()
}
rm(i)

## plot geomean vs. geomean-scaled PC1, for each geo_cluster
for (i in picidae.geo_clusters.to_use) {
  pdf(file=paste("Picidae_community_geomean_vs_geomean_scaled_pgls_PC1_", i, ".pdf", sep=""))
  plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$geomean[grep(i, picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$full_name)], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[grep(i, rownames(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S)),1], ylab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "% of remaining variance)", sep=""), xlab = "Geometric Mean Size", xlim=xylims$geomean.pgls_scaled.pc1.x, ylim=xylims$geomean.pgls_scaled.pc1.y)
  dev.off()
}
rm(i)

## create vectors to use with plotting symbols by species and geo_cluster; note that the numbers are based on the alphabetical order of the My_genus_species values
picidae_species_order <- as.numeric(as.factor(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$My_genus_species))
picidae_geo_cluster_order <- as.numeric(as.factor(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$geo_cluster_shortname))

## plot geomean vs. geomean-scaled PC1, for all populations, with shapes by species and colors by geo_cluster
pdf(file='Picidae_community_geomean_vs_geomean_scaled_PC1_all_geo_colors_taxa_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$geomean, picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1], xlim=xylims$geomean.pgls_scaled.pc1.x, ylim=xylims$geomean.pgls_scaled.pc1.y, xlab = "Geometric mean of measurements", ylab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "% of remaining variance)", sep=""), pch=c(16,0,5,15,6,2,17,1,18,4,3)[picidae_species_order], col=c("blue","forestgreen","cyan1","black","darkgoldenrod","red")[picidae_geo_cluster_order])
dev.off()

## plot geomean vs. geomean-scaled PC1, for all populations, with shapes by geo_clusters
pdf(file='Picidae_community_geomean_vs_geomean_scaled_PC1_all_geo_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$geomean, picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1], xlim=xylims$geomean.pgls_scaled.pc1.x, ylim=xylims$geomean.pgls_scaled.pc1.y, xlab = "Geometric mean of measurements", ylab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "% of remaining variance)", sep=""), pch=c(16,0,5,15,2,17)[picidae_geo_cluster_order])
legend(x="bottomright", pch=c(16,0,5,15,2,17), legend=c("Florida", "Kansas", "Michigan", "Minnesota", "Oregon", "Washington"), cex=0.5)
dev.off()

## plot geomean vs. geomean-scaled PC1, for all populations, with shapes by species
pdf(file='Picidae_community_geomean_vs_geomean_scaled_PC1_all_taxa_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$geomean, picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1], xlim=xylims$geomean.pgls_scaled.pc1.x, ylim=xylims$geomean.pgls_scaled.pc1.y, xlab = "Geometric mean of measurements", ylab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "% of remaining variance)", sep=""), pch=c(16,0,5,15,6,2,17,1,18,4,3)[picidae_species_order])
legend(x="bottomright", pch=c(16,0,5,15,6,2,17,1,18,4,3), legend=gsub("_", " ", gsub("_", " ", unique(sort(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$My_genus_species)))), cex=0.5)
dev.off()

## plot geomean-scaled PC1 vs. geomean-scaled PC2, for all populations, with shapes by species and colors by geo_cluster
pdf(file='Picidae_community_geomean_scaled_PC1_vs_PC2_all_geo_colors_taxa_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,2], xlab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[2], 2), "%)", sep=""), xlim=xylims$pgls_scaled.pc12.x, ylim=xylims$pgls_scaled.pc12.y, pch=c(16,0,5,15,6,2,17,1,18,4,3)[picidae_species_order], col=c("blue","forestgreen","cyan1","black","darkgoldenrod","red")[picidae_geo_cluster_order])
dev.off()

## plot geomean-scaled PC1 vs. geomean-scaled PC2, for all populations, with shapes by geo_clusters
pdf(file='Picidae_community_geomean_scaled_PC1_vs_PC2_all_geo_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,2], xlab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[2], 2), "%)", sep=""), xlim=xylims$pgls_scaled.pc12.x, ylim=xylims$pgls_scaled.pc12.y, pch=c(16,0,5,15,2,17)[picidae_geo_cluster_order])
legend(x="topright", pch=c(16,0,5,15,2,17), legend=c("Florida", "Kansas", "Michigan", "Minnesota", "Oregon", "Washington"), cex=0.5)
dev.off()

## plot geomean-scaled PC1 vs. geomean-scaled PC2, for all populations, with shapes by species
pdf(file='Picidae_community_geomean_scaled_PC1_vs_PC2_all_taxa_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,1], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$pca$S[,2], xlab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$geomean_scaled_pgls.all_inds$var_percents[2], 2), "%)", sep=""), xlim=xylims$pgls_scaled.pc12.x, ylim=xylims$pgls_scaled.pc12.y, pch=c(16,0,5,15,6,2,17,1,18,4,3)[picidae_species_order])
legend(x="topright", pch=c(16,0,5,15,6,2,17,1,18,4,3), legend=gsub("_", " ", unique(sort(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$My_genus_species))), cex=0.5)
dev.off()

## plot unscaled PC1 vs. unscaled PC2, for all populations, with shapes by species and colors by geo_cluster
pdf(file='Picidae_community_unscaled_PC1_vs_PC2_all_geo_colors_taxa_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,1], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,2], xlab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$var_percents[2], 2), "%)", sep=""), xlim=xylims$unscaled.pc12.x, ylim=xylims$unscaled.pc12.y, pch=c(16,0,5,15,6,2,17,1,18,4,3)[picidae_species_order], col=c("blue","forestgreen","cyan1","black","darkgoldenrod","red")[picidae_geo_cluster_order])
dev.off()

## plot unscaled PC1 vs. unscaled PC2, for all populations, with shapes by geo_clusters
pdf(file='Picidae_community_unscaled_PC1_vs_PC2_all_geo_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,1], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,2], xlab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$var_percents[2], 2), "%)", sep=""), xlim=xylims$unscaled.pc12.x, ylim=xylims$unscaled.pc12.y, pch=c(16,0,5,15,2,17)[picidae_geo_cluster_order])
legend(x="bottomright", pch=c(16,0,5,15,2,17), legend=c("Florida", "Kansas", "Michigan", "Minnesota", "Oregon", "Washington"), cex=0.6)
dev.off()

## plot unscaled PC1 vs. unscaled PC2, for all populations, with shapes by species
pdf(file='Picidae_community_unscaled_PC1_vs_PC2_all_taxa_shapes.pdf',height=6, width=6)
plot(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,1], picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$pca$S[,2], xlab = paste("PC1 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$var_percents[1], 1), "%)", sep=""), ylab = paste("PC2 (", round(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.phyl_pca$unscaled.all_inds$var_percents[2], 2), "%)", sep=""), xlim=xylims$unscaled.pc12.x, ylim=xylims$unscaled.pc12.y, pch=c(16,0,5,15,6,2,17,1,18,4,3)[picidae_species_order])
legend(x="bottomright", pch=c(16,0,5,15,6,2,17,1,18,4,3), legend=gsub("_", " ", gsub("_", " ", unique(sort(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean$all_inds$My_genus_species)))), cex=0.6)
dev.off()


### run tests on data and generate null models

geo_clusters.data_variants.names <- read.csv("Picidae_geo_clusters_data_variants_names.csv", header=TRUE, stringsAsFactors=FALSE) # this file has information on the data variants (unscaled, geomean-scaled, geomean) and the objects they use, so that I can loop over them


### calculate metrics for each population, for geomean, unscaled data, and size-scaled shape data

geo_cluster.comm.metrics <- list()  # generate list to store calculated community metrics

## loop over my data_variants object, calculating community metrics for each community for each data variant
for (i in 1:nrow(geo_clusters.data_variants.names)) {
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n\n")
  data.tmp <- get(geo_clusters.data_variants.names$traits_list_orig[i])[[geo_clusters.data_variants.names$df_name_orig[i]]]
  geo_cluster.comm.metrics[[geo_clusters.data_variants.names$name[i]]] <- list()
  for (j in picidae.geo_clusters.to_use) {
    cat("Starting geo cluster", j, "\n")
    geo_cluster.comm.metrics[[geo_clusters.data_variants.names$name[i]]][[j]] <- calcCommMetrics(traits.obs = data.tmp[data.tmp$geo_cluster_shortname==j,], name_col="My_genus_species", trait_cols=4:ncol(data.tmp))
  }
  cat("\n")
}
rm(i,j,data.tmp)


### now generate null models for each community, and generate distributions of summary statistics based on those null models

# generate list to store null distributions of community metrics
geo_cluster.comm.metrics.null_distributions <- list()

## loop over the data variants, generating expected distribution of community metrics with (1) permutation (2) using all the populations of all the species in my measurements
for (i in 1:nrow(geo_clusters.data_variants.names)) {
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n\n")
  data.tmp <- get(geo_clusters.data_variants.names$traits_list_orig[i])[[geo_clusters.data_variants.names$df_name_orig[i]]]  # get the trait data for the current data variant
  geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["permutation"]][["all_focal_pops"]] <- list()
  for (j in picidae.geo_clusters.to_use) {
    cat("Starting geo cluster", j, "\n")
    geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["permutation"]][["all_focal_pops"]][[j]] <- simulateCommMetrics(nsim=1000, method="permutation", ntaxa=sum(data.tmp$geo_cluster_shortname==j), traits.obs=data.tmp, trait_cols=4:ncol(data.tmp), return.sims=FALSE)
  }
  cat("\n")
}
rm(i,j,data.tmp)


### generate expected distributions of community metrics using phylogenetic null model with BM simulation on the tree

## test BM vs. OU models for phylogenetic models

# generate lists to store results from fitting OU and BM models to data
geo_cluster.comm.phylo_models <- list()
geo_cluster.comm.phylo_models[["BM"]] <- list()
geo_cluster.comm.phylo_models[["OU"]] <- list()

## fit BM and OU models for the whole data set (all populations), and one for each community
for (i in geo_clusters.data_variants.names$name) {  # loop over data variants
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n")
  data.tmp <- get(geo_clusters.data_variants.names$traits_list_orthog[i])[[geo_clusters.data_variants.names$df_name_orthog[i]]]  # get data for current data variant
  
  # get trait matrix from current data object (necessary because the structures vary)
  if (is.data.frame(data.tmp)) {
    data.tmp.matrix <- data.matrix(data.tmp[,4:ncol(data.tmp)])
    rownames(data.tmp.matrix) <- data.tmp$full_name
  } else if ((class(data.tmp) == "list") & (class(data.tmp$pca) == "phyl.pca")) {
    data.tmp.matrix <- data.tmp$pca$S
  } else cat("Data format not recognized")
  
  treedata.tmp <- treedata(phy=multi2di(picidae.RAxML.all.BEAST_calibrated.geo_clusters), data=data.tmp.matrix, warnings=FALSE)  # generate treedata object for all populations
  
  # generate lists to store BM and OU models for current data variant
  geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]] <- list()
  geo_cluster.comm.phylo_models[["OU"]][[geo_clusters.data_variants.names$name[i]]] <- list()
  
  # fit BM and OU model to data from all populations
  geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][["all_pops"]] <- fitContinuous(phy=treedata.tmp$phy, dat=treedata.tmp$dat, model="BM")
  geo_cluster.comm.phylo_models[["OU"]][[geo_clusters.data_variants.names$name[i]]]$all_pops <- fitContinuous(phy=treedata.tmp$phy, dat=treedata.tmp$dat, model="OU")
  
  # loop over geo clusters for current data variant
  for (j in picidae.geo_clusters.to_use) {
    cat("Starting geo cluster", j, "\n")
    treedata.comm.tmp <- treedata(phy=treedata.tmp$phy, data=treedata.tmp$data[grep(j, rownames(treedata.tmp$data)),])  # generate treedata object for current geo cluster for current data variant
    
    # fit BM and OU models to current geo cluster for current data variant
    geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][[j]] <- fitContinuous(phy=treedata.comm.tmp$phy, dat=treedata.comm.tmp$data, model="BM")
    geo_cluster.comm.phylo_models[["OU"]][[geo_clusters.data_variants.names$name[i]]][[j]] <- fitContinuous(phy=treedata.comm.tmp$phy, dat=treedata.comm.tmp$data, model="OU")
  }  
}
rm(i,j,data.tmp, data.tmp.matrix, treedata.tmp, treedata.comm.tmp)

# output models for comparison
for (i in names(geo_cluster.comm.phylo_models[["BM"]])) {  # loop over data variants
  cat("Data variant:", i, "\n\n")
  for (j in names(geo_cluster.comm.phylo_models[["BM"]][[i]])) {  # loop over geo clusters
    cat("\nGeo cluster:", j, "\n\n")
    if ("lik" %in% names(geo_cluster.comm.phylo_models[["BM"]][[i]][[j]])) {
      cat("geomean:    BM: ", geo_cluster.comm.phylo_models[["BM"]][[i]][[j]]$opt$aicc, "    OU: ", geo_cluster.comm.phylo_models[["OU"]][[i]][[j]]$opt$aicc, "\n", sep="")  # print AICc for BM and OU
    } else for (k in names(geo_cluster.comm.phylo_models[["BM"]][[i]][[j]])) {  # if multiple axes, loop over axes
      cat(k, ":    BM: ", geo_cluster.comm.phylo_models[["BM"]][[i]][[j]][[k]]$opt$aicc, "    OU: ", geo_cluster.comm.phylo_models[["OU"]][[i]][[j]][[k]]$opt$aicc, "\n", sep="")  # print AICc for BM and OU
    }
  }
  cat("\nDone with data variant", geo_clusters.data_variants.names$name[i], "\n\n\n")
}
rm(i,j,k)


## generate expected distributions using BM simulation
for (i in 1:nrow(geo_clusters.data_variants.names)) {  # loop over data variants
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n\n")
  data.tmp <- get(geo_clusters.data_variants.names$traits_list_orthog[i])[[geo_clusters.data_variants.names$df_name_orthog[i]]]  # get data for current data variant
  
  # get trait matrix from current data object (necessary because the structures vary)
  if (is.data.frame(data.tmp)) {
    data.tmp.matrix <- data.matrix(data.tmp[,4:ncol(data.tmp)])
    rownames(data.tmp.matrix) <- data.tmp$full_name
  } else if ((class(data.tmp) == "list") & (class(data.tmp$pca) == "phyl.pca")) {
    data.tmp.matrix <- data.tmp$pca$S
  } else cat("Data format not recognized")
  
  geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["phylogenetic"]][["comm_members"]] <- list()
  
  for (j in picidae.geo_clusters.to_use) {  # loop over geo clusters
    cat("Starting geo cluster", j, "\n")
    
    # get the BM sigmasquared values (need to loop over the fitContinuous list in order to get each one)
    treedata.tmp <- treedata(phy=picidae.RAxML.all.BEAST_calibrated.geo_clusters, data=data.tmp.matrix[grep(j, rownames(data.tmp.matrix)),])  # generate a treedata object for the current geo cluster and current data variant (this is done primarily to trim the tree to include only populations from the current geo cluster)
    
    # get the sigma^2 value(s) for the current data variant and current geo cluster
    if ("lik" %in% names(geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][[j]])) {
      sigmasq <- geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][[j]]$opt$sigsq
    } else {
      sigmasq <- numeric()
      for (k in 1:length(geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][[j]])) {
        sigmasq[k] <- geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][[j]][[k]]$opt$sigsq
      }
    }
    
    geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["phylogenetic"]][["comm_members"]][[j]] <- simulateCommMetrics(nsim=1000, method="BM", phy=treedata.tmp$phy, sigmasq=sigmasq, return.sims=FALSE)  # generate the expected distributions of community metrics
  }
}
rm(i,j,k, data.tmp, data.tmp.matrix, treedata.tmp, sigmasq)


### generate simulated communities using permutation and simulated trait evolution from all species

## generate simulated communities using permutation samples of all species
for (i in 1:nrow(geo_clusters.data_variants.names)) {  # loop over data variants
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n\n")
  data.tmp <- get(geo_clusters.data_variants.names$species_traits_list_orig[i])[[geo_clusters.data_variants.names$species_df_name_orig[i]]]  # get data for current data variant
  geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["permutation"]][["all_species"]] <- list()
  
  for (j in picidae.geo_clusters.to_use) {  # loop over geo clusters
    cat("Starting geo cluster", j, "\n")
    geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["permutation"]][["all_species"]][[j]] <- simulateCommMetrics(nsim=1000, method="permutation", ntaxa=sum(get(geo_clusters.data_variants.names$traits_list_orig[i])[[geo_clusters.data_variants.names$df_name_orig[i]]]$geo_cluster_shortname==j), traits.obs=data.tmp, trait_cols=2:ncol(data.tmp), return.sims=FALSE)  # generate expected distribution for current geo cluster and data variant
  }  
}
rm(i,j,data.tmp)

## fit BM models to the phylogeny and species means for all species
for (i in 1:nrow(geo_clusters.data_variants.names)) {  # loop over data variants
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n\n")
  data.tmp <- get(geo_clusters.data_variants.names$species_traits_list_orthog[i])[[geo_clusters.data_variants.names$species_df_name_orthog[i]]]  # get data for current data variant
  
  # get trait matrix from current data object (necessary because the structures vary)
  if (is.data.frame(data.tmp)) {
    data.tmp.matrix <- data.matrix(data.tmp[,2:ncol(data.tmp)])
    rownames(data.tmp.matrix) <- data.tmp$My_genus_species
  } else if (class(data.tmp) == "matrix") {
    data.tmp.matrix <- data.tmp
  } else cat("Data format not recognized")
  
  treedata.tmp <- treedata(phy=multi2di(picidae.RAxML.all.BEAST_calibrated.geo_clusters), data=data.tmp.matrix, warnings=FALSE)  # generate treedata object for current data variant
  
  # fit BM and OU models for current data variant
  geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][["all_species"]] <- fitContinuous(phy=treedata.tmp$phy, dat=treedata.tmp$dat, model="BM")
  geo_cluster.comm.phylo_models[["OU"]][[geo_clusters.data_variants.names$name[i]]][["all_species"]] <- fitContinuous(phy=treedata.tmp$phy, dat=treedata.tmp$dat, model="OU")
}
rm(i,data.tmp, data.tmp.matrix, treedata.tmp)

## generate expected distributions of community metrics using BM simulation of the phylogeny of all species and the fit values for sigmasq
# basically, this uses a BM model fit to all the species means, takes a subset of the species (with n = number of taxa in the focal community), and simulates BM on a tree of those taxa
for (i in 1:nrow(geo_clusters.data_variants.names)) {  # loop over data variants
  cat("Starting data variant", geo_clusters.data_variants.names$name[i], "\n\n")
  data.tmp <- get(geo_clusters.data_variants.names$species_traits_list_orthog[i])[[geo_clusters.data_variants.names$species_df_name_orthog[i]]]  # get data for current data variant
  
  # get trait matrix from current data object (necessary because the structures vary)
  if (is.data.frame(data.tmp)) {
    data.tmp.matrix <- data.matrix(data.tmp[,2:ncol(data.tmp)])
    rownames(data.tmp.matrix) <- data.tmp$My_genus_species
  } else if (class(data.tmp) == "matrix") {
    data.tmp.matrix <- data.tmp
  } else cat("Data format not recognized")

  treedata.tmp <- treedata(phy=picidae.RAxML.all.BEAST_calibrated.geo_clusters, data=data.tmp.matrix)  # generate treedata object for current data variant
  
  # get sigma^2 value(s) for current data variant
  if ("lik" %in% names(geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][["all_species"]])) {
    sigmasq <- geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][["all_species"]]$opt$sigsq
  } else {
    sigmasq <- numeric()
    for (k in 1:length(geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][["all_species"]])) {
      sigmasq[k] <- geo_cluster.comm.phylo_models[["BM"]][[geo_clusters.data_variants.names$name[i]]][["all_species"]][[k]]$opt$sigsq
    }
  }
  
  geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["phylogenetic"]][["all_species"]] <- list()
  
  # loop over geo clusters and generate expected distributions of community metrics for each
  for (j in picidae.geo_clusters.to_use) {
    cat("Starting geo cluster", j, "\n")
    geo_cluster.comm.metrics.null_distributions[[geo_clusters.data_variants.names$name[i]]][["phylogenetic"]][["all_species"]][[j]] <- simulateCommMetrics(nsim=1000, method="BM", phy=treedata.tmp$phy, sigmasq=sigmasq, ntaxa=sum(get(geo_clusters.data_variants.names$traits_list_orig[i])[[geo_clusters.data_variants.names$df_name_orig[i]]]$geo_cluster_shortname==j), return.sims=FALSE)
  }
}
rm(i,j,k, data.tmp, data.tmp.matrix, treedata.tmp, sigmasq)


### compare empirical values of community metrics to expected distributions based on simulations
## need to do comparisons for each data variant, for each null model method, for each set of taxa included in the null model, for each community, for each metric
## numbers returned represent the proportion of simulated communities for which the empirical value of the metric EXCEEDS the value for the simulated community

geo_cluster.comm.metrics.percentages <- list()  # generate list to store percentages for all the various combinations

for (i in geo_clusters.data_variants.names$name) {  # loop over data variants
  cat("Starting data variant", i, "\n\n")
  geo_cluster.comm.metrics.percentages[[i]] <- list()
  for (j in names(geo_cluster.comm.metrics.null_distributions[[i]])) {  # loop over simulation method
    geo_cluster.comm.metrics.percentages[[i]][[j]] <- list()
    for (k in names(geo_cluster.comm.metrics.null_distributions[[i]][[j]])) {  # loop over taxa included in simulation method
      geo_cluster.comm.metrics.percentages[[i]][[j]][[k]] <- list()
      for (l in names(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]])) {  # loop over communities
        
        # generate vector to store percentages
        geo_cluster.comm.metrics.percentages[[i]][[j]][[k]][[l]] <- rep(NA, ncol(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]]))
        names(geo_cluster.comm.metrics.percentages[[i]][[j]][[k]][[l]]) <- colnames(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]])
        
        for (m in colnames(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]])) {  # loop over metrics
          empirical.tmp <- geo_cluster.comm.metrics[[i]][[l]][m]
          simulated.tmp <- geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]][,m]
          geo_cluster.comm.metrics.percentages[[i]][[j]][[k]][[l]][m] <- sum(empirical.tmp > simulated.tmp) / length(simulated.tmp)
        }
      }
    }
  }  
}
rm(i,j,k,l,m, simulated.tmp, empirical.tmp)


### select metrics that are most informative
## using pca to do this; metrics that load very similarly on the PCA axes are strongly correlated, and therefore redundant

geo_cluster.comm.metrics.null_distributions.pca_metrics <- list()  # generate list to store pca of metrics

## generate PCA of the different metrics from the simulated communities, for each data variant, for each null model and species pool inclusion combination, for each community
for (i in names(geo_cluster.comm.metrics.null_distributions)) {  # loop over data variants
  geo_cluster.comm.metrics.null_distributions.pca_metrics[[i]] <- list()
  for (j in names(geo_cluster.comm.metrics.null_distributions[[i]])) {  # loop over null model
    geo_cluster.comm.metrics.null_distributions.pca_metrics[[i]][[j]] <- list()
    for (k in names(geo_cluster.comm.metrics.null_distributions[[i]][[j]])) {  # loop over species pool inclusion
      geo_cluster.comm.metrics.null_distributions.pca_metrics[[i]][[j]][[k]] <- list()
      for (l in names(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]])) {  # loop over communities
        
        # get metrics for current combination
        data.tmp <- matrix(NA, nrow=length(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]][[1]]), ncol=length(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]]), dimnames=list(1:length(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]][[1]]), names(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]])))  
        for (m in 1:length(geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]])) {  # loop over metrics
          data.tmp[,m] <- geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]][[m]]
        }
        
        geo_cluster.comm.metrics.null_distributions.pca_metrics[[i]][[j]][[k]][[l]] <- prcomp(data.tmp, center=TRUE, scale.=TRUE)  # calculate principal comoponents for current combination
      }
    }
  }
}
rm(i,j,k,l,m,data.tmp)

## output biplots of the metrics on PCA axes to determine suitable variables to use
for (i in 1:length(geo_cluster.comm.metrics.null_distributions.pca_metrics[[3]][[1]][[1]])) print(ggbiplot(geo_cluster.comm.metrics.null_distributions.pca_metrics[[3]][[1]][[1]][[i]], choices=c(1,3), varname.size=7) + xlim(-3.5,3.5) + ylim(-3.5,3.5))

# store vector of metrics to include
metrics.to_test <- c("FRic", "FEve", "FDiv", "FDis", "min.dist", "mean.nndist", "var.dist", "rat.mst")


### determine variation among null models (permutation vs. phylogenetic and inclusive vs. limited species pools) in the percentages
## test this by conducting ANOVA of the proportions by the null model, species inclusion, and community, to see which one explains the most variation
## first need to transform the proportions using arcsin(sqrt(proportions)) to make distributions approximately normal

## run ANOVA on the various combinations
geo_cluster.comm.metrics.percentages.anova <- list()
for (i in names(geo_cluster.comm.metrics.percentages)) {  # loop over data variants
  geo_cluster.comm.metrics.percentages.anova[[i]] <- list()
  for (j in metrics.to_use) {  # loop over the metrics to be used
    
    # build up a data object to store the data for the current data variant
    if (j %in% names(geo_cluster.comm.metrics.percentages[[i]][[1]][[1]][[1]])) {
      null.model <- character()
      species.pool <- character()
      community <- character()
      prop.transformed <- numeric()
      for (k in names(geo_cluster.comm.metrics.percentages[[i]])) {  # loop over null model type
        for (l in names(geo_cluster.comm.metrics.percentages[[i]][[k]])) {  # loop over species inclusion
          for (m in names(geo_cluster.comm.metrics.percentages[[i]][[k]][[l]])) {  # loop over community
            null.model <- c(null.model, k)
            species.pool <- c(species.pool, l)
            community <- c(community, m)
            prop.transformed <- c(prop.transformed, asin(sqrt(geo_cluster.comm.metrics.percentages[[i]][[k]][[l]][[m]][[j]])))
          }
        }
      }
      
      data.tmp <- data.frame(null.model, species.pool, community, prop.transformed, stringsAsFactors=TRUE)  # combine data for current combination into a single data frame
      geo_cluster.comm.metrics.percentages.anova[[i]][[j]] <- aov(formula=prop.transformed ~ community + null.model/species.pool, data=data.tmp)  # run ANOVA for current combination
    }
  }
}
rm(i,j,k,l,m, data.tmp, null.model, species.pool, community, prop.transformed)

## output results of ANOVA to console
for (j in names(geo_cluster.comm.metrics.percentages.anova[[1]])) {
  for (i in names(geo_cluster.comm.metrics.percentages.anova)) {
    cat(i,j,"\n")
    print(summary(geo_cluster.comm.metrics.percentages.anova[[i]][[j]]))
    cat("\n\n")
  }      
}
rm(i,j)

metrics.to_use <- c("FRic", "FEve", "FDiv", "FDis", "min.dist", "mean.nndist", "rat.mst")  # drop var.dist from the vector of metrics to use, based on results of ANOVA


### export tables and heatmaps of community metric percentages

## generate data frames with community metric percentages
geo_cluster.comm.metrics.percentages.for_table <- list() # create a list to store the data frames for tables
for (i in names(geo_cluster.comm.metrics.percentages)) {  # loop over data variants
  
  # generate a list to store a vector for each metric
  metrics.tmp <- list()
  for (m in metrics.to_use) {
    metrics.tmp[[m]] <- numeric()
  }
  
  # generate vectors to store information about the variants (community, null model)
  community.tmp <- character()
  null_mod_spec_inc.tmp <- character()
  
  for (l in names(geo_cluster.comm.metrics.percentages[[i]][[1]][[1]])) {  # loop over communities
    for (j in names(geo_cluster.comm.metrics.percentages[[i]])) {  # loop over null model
      for (k in names(geo_cluster.comm.metrics.percentages[[i]][[j]])) {  # loop over species inclusion
        community.tmp <- c(community.tmp, picidae.geo_clusters.name_translation[l])  # add current community name to vector of community names
        null_mod_spec_inc.tmp <- c(null_mod_spec_inc.tmp, paste(j,spec_inc.name_translation[k],sep=" with "))  # add current null model to vector of null models
        for (m in metrics.to_use) {  # loop over metrics
          metrics.tmp[[m]] <- c(metrics.tmp[[m]], geo_cluster.comm.metrics.percentages[[i]][[j]][[k]][[l]][[m]])  # get the percentage for the current combination
        }
      }
    }
  }
  
  # create a data frame to store percentages
  geo_cluster.comm.metrics.percentages.for_table[[i]] <- data.frame(community=community.tmp, null_model=null_mod_spec_inc.tmp)  
  
  # loop over metrics, adding vector for that metric to the data frame if there is data for it
  for (m in metrics.to_use) {
    if (length(metrics.tmp[[m]]) > 0) {  # check that metric is included for this data variant
      geo_cluster.comm.metrics.percentages.for_table[[i]] <- cbind(geo_cluster.comm.metrics.percentages.for_table[[i]], metrics.tmp[[m]])
      names(geo_cluster.comm.metrics.percentages.for_table[[i]])[ncol(geo_cluster.comm.metrics.percentages.for_table[[i]])] <- m
    }
  }
}
rm(i,j,k,l,m, community.tmp, null_mod_spec_inc.tmp, metrics.tmp)

## export data frames of percentages to tables for later use
for (i in names(geo_cluster.comm.metrics.percentages.for_table)) {
  write.xlsx(geo_cluster.comm.metrics.percentages.for_table[[i]], file=paste("Picidae_community_metrics_table_", i, ".xlsx", sep=""))
}
rm(i)

## output heat map of each data frame to visualize it
for (i in names(geo_cluster.comm.metrics.percentages.for_table)) {
  pdf(file=paste("Community_metrics_heatmap_", i, ".pdf", sep=""), width=12, height=8)
  par(omi=c(2,0,0,3))
  heatmap.2(data.matrix(geo_cluster.comm.metrics.percentages.for_table[[i]][,3:ncol(geo_cluster.comm.metrics.percentages.for_table[[i]])]), Rowv = FALSE, Colv = FALSE, labRow=paste(geo_cluster.comm.metrics.percentages.for_table[[i]]$community, geo_cluster.comm.metrics.percentages.for_table[[i]]$null_model, sep=" "), key=FALSE, trace="none", col=gray(seq(0,1,by=0.01)))
  dev.off()
}

## output percentages for each metric for each data variant, null model, and species inclusion variant, for all communities together
# sink(file="Picidae_community_trait_distributions_metrics_community_comparison_all.txt")  # use this if creating a single file
for (i in names(geo_cluster.comm.metrics.percentages)) {  # loop over data variants
  sink(file=paste("Picidae_community_trait_distributions_metrics_community_comparison_", i, ".txt", sep=""))
  # cat("Data variant:", i, "\n\n")
  for (j in intersect(names(geo_cluster.comm.metrics.percentages[[i]][[1]][[1]][[1]]), metrics.to_use)) {  # loop over metric
    cat("Metric:", j, "\n\n")
    for (k in names(geo_cluster.comm.metrics.percentages[[i]])) {  # loop over null model type
      for (l in names(geo_cluster.comm.metrics.percentages[[i]][[k]])) {  # loop over species inclusion)
        cat("Null model ", k, "; Species pool ", l, "\n", sep="")
        for (m in names(geo_cluster.comm.metrics.percentages[[i]][[k]][[l]])) {  # loop over communities
          cat(geo_cluster.comm.metrics.percentages[[i]][[k]][[l]][[m]][[j]], "\t", m, "\n", sep="")  # output percentage for current combination
        }
        cat("\n")
      }
    }
    cat("\n")
  }
  cat("\n")
  sink()  # use this if creating separate file for each data variant
}
# sink()  # use this if creating a single file
rm(i,j,k,l,m)


### generate expected distributions using the within-species population permutation (WSPP), and compare these distribution to the empirical values and distributions from the general permutation
## the WSPP uses only the species in the given community, and randomly selects a single population from each of those species

## generate expected distributions using WSPP for two communities of interest (Washington, Florida)
geo_cluster.comm.metrics.wspp.null_distributions <- list()
for (i in names(geo_cluster.comm.metrics)) {  # loop over data variants
  cat("Starting data variant", i, "\n\n")
  geo_cluster.comm.metrics.wspp.null_distributions[[i]] <- list()
  data.tmp <- get(geo_clusters.data_variants.names$traits_list_orig[which(geo_clusters.data_variants.names$name==i)])[[geo_clusters.data_variants.names$df_name_orig[which(geo_clusters.data_variants.names$name==i)]]]  # get data for current data variant
  
  for (j in c("FL_Alachua", "WA_King")) {  # loop over the two communities of interest
    cat("Starting geo cluster", j, "\n")
    data.tmp.comm <- data.tmp[data.tmp$My_genus_species %in% (data.tmp$My_genus_species[data.tmp$geo_cluster_shortname==j]),]  # trim data frame to include only populations of species found in current geo cluster
    geo_cluster.comm.metrics.wspp.null_distributions[[i]][[j]] <- simulateCommMetrics(nsim=1000, method="permutation", ntaxa=sum(data.tmp$geo_cluster_shortname==j), traits.obs=data.tmp.comm, trait_cols=4:ncol(data.tmp.comm), metric=metrics.to_use, return.sims=FALSE)  # generate distribution with WSPP, for the metrics to be used
  }
  cat("\n")
}
rm(i, j, data.tmp, data.tmp.comm)

## get the percentages of the empirical community metrics vs. WSPP (to test if the empirical communities are unusual among possible communities with those species)
geo_cluster.comm.metrics.wspp.percentages <- list()
for (i in names(geo_cluster.comm.metrics.wspp.null_distributions)) {  # loop over data variants
  geo_cluster.comm.metrics.wspp.percentages[[i]] <- list()
  for (j in names(geo_cluster.comm.metrics.wspp.null_distributions[[i]])) {  # loop over communities
    geo_cluster.comm.metrics.wspp.percentages[[i]][[j]] <- numeric()
    for (k in colnames(geo_cluster.comm.metrics.wspp.null_distributions[[i]][[j]])) {  # loop over metrics
      empirical.tmp <- geo_cluster.comm.metrics[[i]][[j]][k]  # get empirical data for current combination
      simulated.tmp <- geo_cluster.comm.metrics.wspp.null_distributions[[i]][[j]][,k]  # get expected distribution from WSPP for current combination
      geo_cluster.comm.metrics.wspp.percentages[[i]][[j]][[k]] <- sum(empirical.tmp > simulated.tmp) / length(simulated.tmp)  # calculate percentage for current combination
    } 
  }
}
rm(i,j,k,empirical.tmp,simulated.tmp)

## generate data frame of percentages for empirical community metrics vs. WSPP
data_variant.tmp <- character()  # initialize vector to store data variant names
community.tmp <- character()  # initialize vector to store community names
metrics.tmp <- list()  # initialize vector to store metric names
for (m in metrics.to_use) {  # generate a set of numeric vectors to store the metrics
  metrics.tmp[[m]] <- numeric()
}
for (i in names(geo_cluster.comm.metrics.wspp.percentages)) {  # loop over data variants
  for (j in names(geo_cluster.comm.metrics.wspp.percentages[[i]])) {  # loop over communities
    data_variant.tmp <- c(data_variant.tmp, i)
    community.tmp <- c(community.tmp, picidae.geo_clusters.name_translation[j])
    for (m in metrics.to_use) {  # loop over metrics
      metrics.tmp[[m]] <- c(metrics.tmp[[m]], geo_cluster.comm.metrics.wspp.percentages[[i]][[j]][m])  # add the current value to the vector of metrics
    }
  }
}
# create the data frame, and add the vectors for each metric to the data frame as additional columns
geo_cluster.comm.metrics.wspp.percentages.for_table <- data.frame(data_variant=data_variant.tmp, community=community.tmp)
for (m in metrics.to_use) {
  if (length(metrics.tmp[[m]]) > 0) {  # check that metric is included for this data variant
    geo_cluster.comm.metrics.wspp.percentages.for_table <- cbind(geo_cluster.comm.metrics.wspp.percentages.for_table, metrics.tmp[[m]])
    names(geo_cluster.comm.metrics.wspp.percentages.for_table)[ncol(geo_cluster.comm.metrics.wspp.percentages.for_table)] <- m
  }
}
rm(i,j,m, community.tmp, metrics.tmp, data_variant.tmp)

## export data frame of percentages for empirical community metrics vs. WSPP as a table for later use
write.xlsx(geo_cluster.comm.metrics.wspp.percentages.for_table, file="Picidae_community_metrics_empirical_vs_WSPP_table_all_data_variants.xlsx")

## output heatmap of percentages for empirical community metrics vs. WSPP
pdf(file="Picidae_community_metrics_empirical_vs_WSPP_heatmap_all_data_variants.pdf", height=5, width=12)
par(omi=c(2,0,0,3))
heatmap.2(data.matrix(geo_cluster.comm.metrics.wspp.percentages.for_table[c(2,4,6,5),3:ncol(geo_cluster.comm.metrics.wspp.percentages.for_table)]), Rowv = FALSE, Colv = FALSE, key=FALSE, trace="none", col=gray(seq(0,1,by=0.01)))
dev.off()

## calculate the percentages of each iteration of the WSPP vs. other null models (general permutation and phylogenetic)
geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions <- list()
for (i in names(geo_cluster.comm.metrics.wspp.null_distributions)) {  # loop over data variants
  geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]] <- list()
  for (j in names(geo_cluster.comm.metrics.null_distributions[[i]])) {  # loop over null model
    geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]] <- list()
    for (k in names(geo_cluster.comm.metrics.null_distributions[[i]][[j]])) {  # loop over species inclusion
      geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]] <- list()
      for (l in names(geo_cluster.comm.metrics.wspp.null_distributions[[i]])) {  # loop over communities
        geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]][[l]] <- list()
        for (m in colnames(geo_cluster.comm.metrics.wspp.null_distributions[[i]][[l]])) {  # loop over metrics
          geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]][[l]][[m]] <- numeric()
          simulated.tmp <- geo_cluster.comm.metrics.null_distributions[[i]][[j]][[k]][[l]][,m]  # get the expected distribution for the current metric using the current general null model
            for (n in 1:nrow(geo_cluster.comm.metrics.wspp.null_distributions[[i]][[l]])) {  # loop over simulations
              empirical.tmp <- geo_cluster.comm.metrics.wspp.null_distributions[[i]][[l]][n,m]  # get the WSPP metric value for the current combination and current iteration
              geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]][[l]][[m]][n] <- sum(empirical.tmp > simulated.tmp) / length(simulated.tmp)  # calculate the percentage for the current iteration
          }
        }
      }
    }
  }
}
rm(i,j,k,l,m,n,empirical.tmp,simulated.tmp)

## calculate median of percentages of the WSPP vs. other null models (general permutation and phylogenetic), and output histograms of the distributions
geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median <- list()
for (i in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions)) {  # loop over data variants
  cat("Starting data variant", i, "\n\n")
  geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]] <- list()
  for (j in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]])) {  # loop over null model
    cat("Starting null model", j, "\n\n")
    geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[j]] <- list()
    
    # output histograms of distribution for current combination
    pdf(file=paste("Picidae_WSPP_vs_other_null_distributions", i, "_", j, "_hists.pdf", sep=""), width=(2*length(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]])*length(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[1]])), height=2*length(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[1]][[1]]))
    par(mfcol=c(length(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[1]][[1]]),length(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]])*length(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[1]])))
    for (k in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]])) {  # loop over species inclusion
      cat("Starting species inclusion", k, "\n\n")
      geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[j]][[k]] <- list()
      for (l in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]])) {  # loop over community
        cat("Starting community", l, "\n")
        geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[j]][[k]][[l]] <- numeric()
        for (m in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]][[l]])) {  # loop over metric
          hist(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]][[l]][[m]], ylab=NULL, xlab=NULL, main=NULL, xlim=c(0,1), breaks=10)  # output histogram
          geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[j]][[k]][[l]][[m]] <- median(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions[[i]][[j]][[k]][[l]][[m]])  # calculate median percentage of the WSPP vs. other null models (general permutation and phylogenetic)
        }
        cat("\n")
      }
    }
    dev.off()
  }  
}
rm(i,j,k,l,m)

## generate data frame with the median percentage of the WSPP community metrics vs. other null distributions
geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table <- list() # create a list to store the data frames for tables for all 3 data variants
for (i in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median)) {  # loop over data variants
  metrics.tmp <- list()
  for (m in metrics.to_use) {  # generate a set of numeric vectors to store the metrics
    metrics.tmp[[m]] <- numeric()
  }
  community.tmp <- character()
  null_mod_spec_inc.tmp <- character()
  for (l in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[1]][[1]])) {  # loop over communities
    for (j in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]])) {  # loop over null model
      for (k in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[j]])) {  # loop over species inclusion
        community.tmp <- c(community.tmp, picidae.geo_clusters.name_translation[l])
        null_mod_spec_inc.tmp <- c(null_mod_spec_inc.tmp, paste(j,spec_inc.name_translation[k],sep=" with "))
        for (m in metrics.to_use) {  # loop over metrics
          metrics.tmp[[m]] <- c(metrics.tmp[[m]], geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median[[i]][[j]][[k]][[l]][m])
        }
      }
    }
  }
  geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]] <- data.frame(community=community.tmp, null_model=null_mod_spec_inc.tmp)
  
  # add columns to data frame for each metric included
  for (m in metrics.to_use) {
    if (length(metrics.tmp[[m]]) > 0) {  # check that metric is included for this data variant
      geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]] <- cbind(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]], metrics.tmp[[m]])
      names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]])[ncol(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]])] <- m
    }
  }
}
rm(i,j,k,l,m, community.tmp, null_mod_spec_inc.tmp, metrics.tmp)

## export table of the median percentage of the WSPP community metrics vs. other null distributions
for (i in names(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table)) {
  write.xlsx(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]], file=paste("Picidae_WSPP_vs_other_null_distributions_medians_table_", i, ".xlsx", sep=""))
}
rm(i)

## output heatmaps of the median percentage of the WSPP community metrics vs. other null distributions
# generate matrix to make heatmap
geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap <- geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[1]]
for (i in 2:3) {
  geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap <- rbind(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap, geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_table[[i]])
}
# output heatmap
pdf(file="Picidae_WSPP_vs_other_null_distributions_medians_heatmap.pdf", width=12, height=7)
par(omi=c(2,0,0,3))
heatmap.2(data.matrix(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap[c(5:8,13:16,21:24,17:20),3:ncol(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap)]), Rowv = FALSE, Colv = FALSE, labRow=paste(geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap$community[c(5:8,13:16,21:24,17:20)], geo_cluster.comm.metrics.wspp.percentages_of_other_null_distributions.median.for_heatmap$null_model[c(5:8,13:16,21:24,17:20)], sep=" "), key=FALSE, trace="none", col=gray(seq(0,1,by=0.01)))
dev.off()


### test displacement towards/away from nearest neighbors in trait space

## calculate the geomean for the mean of all populations of each species
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops.geomean <- list()
for (i in names(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops)) {
  picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops.geomean[[i]] <- data.frame(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops[[i]][1], geomean=apply(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops[[i]][,2:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops[[i]])],MARGIN=1,mean))
}
rm(i)

## calculate the residuals from PGLS regression on geomean for the mean of all populations of each species (using the PGLS models generated by from the population means)
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops.geomean_scaled <- list()
picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops.geomean_scaled$all_inds.pgls <- sizeScale.pgls.noPhylo(model.pgls.traits_geomean=picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.geomean_scaled$all_inds.pgls_fits$model.pgls.traits_geomean, data.new=picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops$all_inds, scale_cols=2:ncol(picidae.morph_combined.log.fully_reduced.imputed.geo_clusters.means_of_all_pops$all_inds))


### calculate the displacement of each population in each community

geo_cluster.pops.displacement <- list()
for (i in geo_clusters.data_variants.names$name) {  # loop over data variants
  cat("\nStarting data variant", i, "\n\n")
  geo_cluster.pops.displacement[[i]] <- list()
  traits.obs.tmp <- get(geo_clusters.data_variants.names$traits_list_orig[which(geo_clusters.data_variants.names$name==i)])[[geo_clusters.data_variants.names$df_name_orig[which(geo_clusters.data_variants.names$name==i)]]]
  traits.tax_avg.tmp <- get(geo_clusters.data_variants.names$pop_means_traits_list_orig[which(geo_clusters.data_variants.names$name==i)])[[geo_clusters.data_variants.names$pop_means_df_name_orig[which(geo_clusters.data_variants.names$name==i)]]]
  for (k in c("nndist", "gravity")) {  # loop over versions of displacement (nndist and gravity)
    geo_cluster.pops.displacement[[i]][[k]] <- list()
    for (j in picidae.geo_clusters.to_use) {  # loop over communities
      cat("Starting geo cluster", j, "\n")
      traits.obs.tmp.comm <- data.matrix(traits.obs.tmp[traits.obs.tmp$geo_cluster_shortname==j,4:ncol(traits.obs.tmp)])
      rownames(traits.obs.tmp.comm) <- traits.obs.tmp$My_genus_species[traits.obs.tmp$geo_cluster_shortname==j]
      traits.tax_avg.tmp.comm <- data.matrix(traits.tax_avg.tmp[traits.tax_avg.tmp$My_genus_species %in% rownames(traits.obs.tmp.comm), 2:ncol(traits.tax_avg.tmp)])
      rownames(traits.tax_avg.tmp.comm) <- traits.tax_avg.tmp$My_genus_species[traits.tax_avg.tmp$My_genus_species %in% rownames(traits.obs.tmp.comm)]
      if (k=="nndist") {
        geo_cluster.pops.displacement[[i]][[k]][[j]] <- calcDisplacementNNDist(traits.obs = traits.obs.tmp.comm, traits.tax_avg=traits.tax_avg.tmp.comm)
      } else if (k=="gravity") {
        geo_cluster.pops.displacement[[i]][[k]][[j]] <- calcDisplacementGravity(traits.obs = traits.obs.tmp.comm, traits.tax_avg=traits.tax_avg.tmp.comm)
      }
    }
  }
}
rm(i,j, k, traits.obs.tmp, traits.tax_avg.tmp, traits.obs.tmp.comm, traits.tax_avg.tmp.comm)

## test if the displacements by community have more negative or positive values than expected under a binomial distribution with p=0.5
geo_cluster.pops.displacement.comm.pbinom <- list()
for (i in names(geo_cluster.pops.displacement)) {
  geo_cluster.pops.displacement.comm.pbinom[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {
    for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {
      data.tmp <- na.omit(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement)
      geo_cluster.pops.displacement.comm.pbinom[[i]][[k]][[j]] <- list()
      geo_cluster.pops.displacement.comm.pbinom[[i]][[k]][[j]]$negative <- pbinom(q=sum(data.tmp>0), size=length(data.tmp), prob=0.5)
      geo_cluster.pops.displacement.comm.pbinom[[i]][[k]][[j]]$positive <- pbinom(q=sum(data.tmp<0), size=length(data.tmp), prob=0.5)
    }
  }
}
rm(i,j,k, data.tmp)

# test if the displacements by species have more negative or positive values than expected under a binomial distribution with p=0.5
geo_cluster.pops.displacement.spec.pbinom <- list()
for (i in names(geo_cluster.pops.displacement)) {  # loop over data variant
  cat("Starting data variant", i, "\n")
  geo_cluster.pops.displacement.spec.pbinom[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {  # loop over displacement method
    cat("Starting method", k, "\n")
    geo_cluster.pops.displacement.spec.pbinom[[i]][[k]] <- list()
    for (l in unique(get(geo_clusters.data_variants.names$traits_list_orig[which(geo_clusters.data_variants.names$name==i)])[[geo_clusters.data_variants.names$df_name_orig[which(geo_clusters.data_variants.names$name==i)]]]$My_genus_species)) {  # loop over species
      data.tmp <- numeric()
      for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {  # loop over communities
        data.tmp <- c(data.tmp, na.omit(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement[geo_cluster.pops.displacement[[i]][[k]][[j]]$taxon.fp==l]))
      }
      if (length(data.tmp) > 1) {
        geo_cluster.pops.displacement.spec.pbinom[[i]][[k]][[l]][["negative"]] <- pbinom(q=sum(data.tmp>0), size=length(data.tmp), prob=0.5)
        geo_cluster.pops.displacement.spec.pbinom[[i]][[k]][[l]][["positive"]] <- pbinom(q=sum(data.tmp<0), size=length(data.tmp), prob=0.5)
      } else {
        geo_cluster.pops.displacement.spec.pbinom[[i]][[k]][[l]][["negative"]] <- NA
        geo_cluster.pops.displacement.spec.pbinom[[i]][[k]][[l]][["positive"]] <- NA
      }
    }
    cat("\n")
  }
}
rm(i,j,k,l,data.tmp)

# test total number of negatives and total number of positives across all communities against binomical distribution
geo_cluster.pops.all_displacements.pbinom <- list()
for (i in names(geo_cluster.pops.displacement)) {
  geo_cluster.pops.all_displacements.pbinom[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {
    total.tmp <- 0
    negative.tmp <- 0
    positive.tmp <- 0
    for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {
      data.tmp <- na.omit(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement)
      total.tmp <- total.tmp + length(data.tmp)
      negative.tmp <- negative.tmp + sum(data.tmp<0)
      positive.tmp <- positive.tmp + sum(data.tmp>0)
    }
    geo_cluster.pops.all_displacements.pbinom[[i]][[k]][["negative"]] <- pbinom(q=positive.tmp, size=total.tmp, prob=0.5)
    geo_cluster.pops.all_displacements.pbinom[[i]][[k]][["positive"]] <- pbinom(q=negative.tmp, size=total.tmp, prob=0.5)
  }
}
rm(i,j,k, data.tmp, total.tmp, negative.tmp, positive.tmp)

## function to calculate standard error
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(na.omit(x)))

## calculate the mean displacement of the populations in a given community
geo_cluster.pops.displacement.comm.mean <- list()
for (i in names(geo_cluster.pops.displacement)) {
  geo_cluster.pops.displacement.comm.mean[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {  # loop over displacement method
    geo_cluster.pops.displacement.comm.mean[[i]][[k]] <- list()
    geo_cluster.pops.displacement.comm.mean[[i]][[k]][["mean"]] <- numeric()
    geo_cluster.pops.displacement.comm.mean[[i]][[k]][["se"]] <- numeric()
    for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {
      geo_cluster.pops.displacement.comm.mean[[i]][[k]][["mean"]][[j]] <- mean(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement, na.rm=TRUE)
      geo_cluster.pops.displacement.comm.mean[[i]][[k]][["se"]][[j]] <- se(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement)
    }
  }
}
rm(i,j,k)

# test if these mean displacements by community are significantly different from 0, assuming normally distributed values (need to update this to use pt and se rather than pnorm and sd)
geo_cluster.pops.displacement.comm.mean.t_test <- list()
for (i in names(geo_cluster.pops.displacement)) {
  geo_cluster.pops.displacement.comm.mean.t_test[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {
    for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {
      geo_cluster.pops.displacement.comm.mean.t_test[[i]][[k]][[j]] <- t.test(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement, na.rm=TRUE)
    }
  }
}
rm(i,j,k)

## calculate the mean displacement of the populations of a given species, and conduct t-tests on those means
geo_cluster.pops.displacement.spec.mean <- list()
geo_cluster.pops.displacement.spec.mean.t_test <- list()
for (i in names(geo_cluster.pops.displacement)) {  # loop over data variant
  cat("Starting data variant", i, "\n")
  geo_cluster.pops.displacement.spec.mean[[i]] <- list()
  geo_cluster.pops.displacement.spec.mean.t_test[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {  # loop over displacement method
    cat("Starting method", k, "\n")
    geo_cluster.pops.displacement.spec.mean[[i]][[k]] <- list()
    geo_cluster.pops.displacement.spec.mean[[i]][[k]][["mean"]] <- numeric()
    geo_cluster.pops.displacement.spec.mean[[i]][[k]][["se"]] <- numeric()
    geo_cluster.pops.displacement.spec.mean.t_test[[i]][[k]] <- list()
    for (l in unique(get(geo_clusters.data_variants.names$traits_list_orig[which(geo_clusters.data_variants.names$name==i)])[[geo_clusters.data_variants.names$df_name_orig[which(geo_clusters.data_variants.names$name==i)]]]$My_genus_species)) {  # loop over species
      data.tmp <- numeric()
      for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {  # loop over communities
        data.tmp <- c(data.tmp, na.omit(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement[geo_cluster.pops.displacement[[i]][[k]][[j]]$taxon.fp==l]))
      }
      if (length(data.tmp >= 2)) {
      geo_cluster.pops.displacement.spec.mean[[i]][[k]][["mean"]][[l]] <- mean(data.tmp, na.rm=TRUE)
      geo_cluster.pops.displacement.spec.mean[[i]][[k]][["se"]][[l]] <- se(data.tmp)
      geo_cluster.pops.displacement.spec.mean.t_test[[i]][[k]][[l]] <- t.test(data.tmp, na.rm=TRUE)
      } else {
        geo_cluster.pops.displacement.spec.mean[[i]][[k]][["mean"]][[l]] <- NA
        geo_cluster.pops.displacement.spec.mean[[i]][[k]][["se"]][[l]] <- NA
        geo_cluster.pops.displacement.spec.mean.t_test[[i]][[k]][[l]] <- NA
      }
    }
    cat("\n")
  }
}
rm(i,j,k,l,data.tmp)

## get the mean and t-test probability across all populations for each data variant
geo_cluster.pops.all_displacements.mean <- list()
geo_cluster.pops.all_displacements.mean.t_test <- list()
for (i in names(geo_cluster.pops.displacement)) {  # loop over data variant
  geo_cluster.pops.all_displacements.mean[[i]] <- list()
  geo_cluster.pops.all_displacements.mean.t_test[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement[[i]])) {  # loop over method for calculating displacement
    data.tmp <- numeric()
    geo_cluster.pops.all_displacements.mean[[i]][[k]] <- numeric()
    for (j in names(geo_cluster.pops.displacement[[i]][[k]])) {  # loop over communities
      data.tmp <- c(data.tmp, na.omit(geo_cluster.pops.displacement[[i]][[k]][[j]]$displacement))
    }
    geo_cluster.pops.all_displacements.mean[[i]][[k]][["mean"]] <- mean(data.tmp, na.rm=TRUE)
    geo_cluster.pops.all_displacements.mean[[i]][[k]][["se"]] <- se(data.tmp)
    geo_cluster.pops.all_displacements.mean.t_test[[i]][[k]] <- t.test(data.tmp, na.rm=TRUE)
  }
}
rm(i,j,k, data.tmp)

## this generate null models for displacement by building communities made up of a random selection of the populations from the species in the given community
geo_cluster.pops.displacement.simulations <- list()
for (i in 1:length(geo_clusters.data_variants.names$name)) {
  cat("\nStarting data variant", geo_clusters.data_variants.names$name[i], "\n")
  geo_cluster.pops.displacement.simulations[[geo_clusters.data_variants.names$name[i]]] <- list()
  traits.obs.tmp <- get(geo_clusters.data_variants.names$traits_list_orig[i])[[geo_clusters.data_variants.names$df_name_orig[i]]]
  traits.tax_avg.tmp <- get(geo_clusters.data_variants.names$pop_means_traits_list_orig[i])[[geo_clusters.data_variants.names$pop_means_df_name_orig[i]]]
  geo_cluster.pops.displacement.simulations[[geo_clusters.data_variants.names$name[i]]][["gravity"]] <- geo_cluster.pops.displacement.simulations[[geo_clusters.data_variants.names$name[i]]][["nndist"]] <- list()
  for (j in picidae.geo_clusters.to_use) {
    cat("\n\nStarting geo cluster", j, "\n")
    geo_cluster.pops.displacement.simulations[[geo_clusters.data_variants.names$name[i]]][["nndist"]][[j]] <- simulate.displacement(nsim=1000, taxa=traits.obs.tmp$My_genus_species[traits.obs.tmp$geo_cluster_shortname==j], traits.obs=traits.obs.tmp, traits.tax_avg=traits.tax_avg.tmp, method="nndist")
    geo_cluster.pops.displacement.simulations[[geo_clusters.data_variants.names$name[i]]][["gravity"]][[j]] <- simulate.displacement(nsim=1000, taxa=traits.obs.tmp$My_genus_species[traits.obs.tmp$geo_cluster_shortname==j], traits.obs=traits.obs.tmp, traits.tax_avg=traits.tax_avg.tmp, method="gravity")
  }
}
rm(i,j,traits.obs.tmp,traits.tax_avg.tmp)

## calculate the mean values of displacement from each of the simulated communities, to generate a distribution of mean values for each data variant and community
geo_cluster.pops.displacement.mean.distributions <- list()
for (i in names(geo_cluster.pops.displacement.simulations)) {
  cat("\nStarting data variant", i, "\n\n")
  geo_cluster.pops.displacement.mean.distributions[[i]] <- list()
  for (l in names(geo_cluster.pops.displacement.simulations[[i]])) {
    geo_cluster.pops.displacement.mean.distributions[[i]][[l]] <- list()
    for (j in names(geo_cluster.pops.displacement.simulations[[i]][[l]])) {
      cat("\n\nStarting geo cluster", j, "\n")
      geo_cluster.pops.displacement.mean.distributions[[i]][[l]][[j]] <- rep(NA, (length(geo_cluster.pops.displacement.simulations[[i]][[l]][[j]])))
      for (k in 1:length(geo_cluster.pops.displacement.simulations[[i]][[l]][[j]])) {
        geo_cluster.pops.displacement.mean.distributions[[i]][[l]][[j]][k] <- mean(geo_cluster.pops.displacement.simulations[[i]][[l]][[j]][[k]]$displacement, na.rm=TRUE)
      }
    }
  }
}
rm(i,j,k,l)

## get the significance value (as the % of simulations that have mean displacements greater than the observed mean displacement)
geo_cluster.pops.displacement.mean.percentages <- list()
for (i in names(geo_cluster.pops.displacement.mean)) {
  cat("\nStarting data variant", i, "\n\n")
  geo_cluster.pops.displacement.mean.percentages[[i]] <- list()
  for (k in names(geo_cluster.pops.displacement.mean[[i]])) {
    for (j in names(geo_cluster.pops.displacement.mean[[i]][[k]])) {
      cat("\n\nStarting geo cluster", j, "\n")
      geo_cluster.pops.displacement.mean.percentages[[i]][[k]][[j]] <- sum(geo_cluster.pops.displacement.mean[[i]][[k]][[j]] > geo_cluster.pops.displacement.mean.distributions[[i]][[k]][[j]]) / length(geo_cluster.pops.displacement.mean.distributions[[i]][[k]][[j]])
    }
  }
}
rm(i,j,k)
